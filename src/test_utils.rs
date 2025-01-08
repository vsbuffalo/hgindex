// test_utils.rs

#[cfg(test)]
pub mod test_utils {
    use std::env;
    use std::path::{Path, PathBuf};
    use tempfile;

    pub struct TestDir {
        dir: PathBuf,
        #[allow(dead_code)] // CHECK: this seems like warning is in error
        temp_dir: Option<tempfile::TempDir>,
    }

    impl TestDir {
        pub fn new(prefix: &str) -> std::io::Result<Self> {
            let keep_output = env::var("KEEP_TEST_OUTPUT").is_ok();
            if keep_output {
                let output_dir = env::current_dir()?.join("test_output").join(prefix);
                std::fs::create_dir_all(&output_dir)?;
                Ok(TestDir {
                    dir: output_dir,
                    temp_dir: None,
                })
            } else {
                let temp_dir = tempfile::tempdir()?;
                let dir = temp_dir.path().to_path_buf();
                Ok(TestDir {
                    dir,
                    temp_dir: Some(temp_dir),
                })
            }
        }

        pub fn path(&self) -> &Path {
            &self.dir
        }
    }
}
