// test_utils.rs

#[cfg(test)]
pub mod test_utils {
    use std::env;
    use std::path::{Path, PathBuf};

    pub struct TestDir {
        dir: PathBuf,
        _temp_dir: Option<tempfile::TempDir>, // Some if temporary, None if persistent
    }

    impl TestDir {
        pub fn new(prefix: &str) -> std::io::Result<Self> {
            // Check if KEEP_TEST_OUTPUT environment variable is set
            let keep_output = env::var("KEEP_TEST_OUTPUT").is_ok();

            if keep_output {
                // Create a directory in the current directory or target/debug
                let output_dir = env::current_dir()?.join("test_output").join(prefix);
                std::fs::create_dir_all(&output_dir)?;
                Ok(TestDir {
                    dir: output_dir,
                    _temp_dir: None,
                })
            } else {
                // Use tempfile as before
                let temp_dir = tempfile::tempdir()?;
                let dir = temp_dir.path().to_path_buf();
                Ok(TestDir {
                    dir,
                    _temp_dir: Some(temp_dir),
                })
            }
        }

        pub fn path(&self) -> &Path {
            &self.dir
        }
    }
}
