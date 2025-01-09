// io.rs

#[cfg(feature = "cli")]
pub mod io {
    use flate2::read::GzDecoder;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::fs::File;
    use std::io::{self, BufRead, BufReader, BufWriter, Error, Read, Seek, Write};
    use std::path::{Path, PathBuf};
    use thiserror::Error;

    const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];
    const DEFAULT_BUFFER_SIZE: usize = 128 * 1024;

    #[derive(Error, Debug)]
    pub enum IoError {
        #[error("IO error: {0}")]
        IoError(#[from] std::io::Error),
        #[error("Invalid or corrupted gzip header")]
        InvalidGzipHeader,
    }

    pub struct InputStream {
        filepath: PathBuf,
    }

    impl InputStream {
        pub fn new(filepath: &Path) -> Self {
            Self {
                filepath: filepath.into(),
            }
        }

        fn is_gzipped(file: &mut File) -> Result<bool, IoError> {
            let mut header = [0u8; 2];
            // Read the first two bytes
            file.read_exact(&mut header)?;
            // Reset the file pointer
            file.rewind()?;
            Ok(header == GZIP_MAGIC)
        }

        pub fn reader(&self) -> Result<BufReader<Box<dyn Read>>, IoError> {
            let mut file = File::open(&self.filepath)?;
            let reader: Box<dyn Read> = if Self::is_gzipped(&mut file)? {
                Box::new(GzDecoder::new(file))
            } else {
                Box::new(file)
            };
            Ok(BufReader::with_capacity(DEFAULT_BUFFER_SIZE, reader))
        }

        pub fn has_header(&self, expect: &str) -> Result<bool, IoError> {
            let mut buf_reader = self.reader()?;
            let mut first_line = String::new();
            buf_reader.read_line(&mut first_line)?;
            Ok(first_line.starts_with(expect))
        }
    }

    #[derive(Clone)]
    pub struct OutputStreamBuilder {
        filepath: Option<PathBuf>,
        buffer_size: usize,
        compression_level: Compression,
    }

    impl Default for OutputStreamBuilder {
        fn default() -> Self {
            Self {
                filepath: None,
                buffer_size: DEFAULT_BUFFER_SIZE,
                compression_level: Compression::default(),
            }
        }
    }

    impl OutputStreamBuilder {
        pub fn new() -> Self {
            Self::default()
        }

        pub fn filepath(mut self, path: Option<impl AsRef<Path>>) -> Self {
            self.filepath = path.map(|p| p.as_ref().to_path_buf());
            self
        }

        pub fn buffer_size(mut self, size: usize) -> Self {
            self.buffer_size = size;
            self
        }

        pub fn compression_level(mut self, level: Option<Compression>) -> Self {
            if let Some(level) = level {
                self.compression_level = level;
            } else {
                self.compression_level = Compression::best();
            }
            self
        }

        pub fn build(self) -> OutputStream {
            OutputStream {
                filepath: self.filepath,
                buffer_size: self.buffer_size,
                compression_level: self.compression_level,
            }
        }
    }

    pub struct OutputStream {
        filepath: Option<PathBuf>,
        buffer_size: usize,
        compression_level: Compression,
    }

    impl OutputStream {
        pub fn new(filepath: Option<impl AsRef<Path>>) -> Self {
            OutputStreamBuilder::new().filepath(filepath).build()
        }

        pub fn builder() -> OutputStreamBuilder {
            OutputStreamBuilder::new()
        }

        fn should_compress(&self) -> bool {
            self.filepath
                .as_ref()
                .map_or(false, |p| p.extension().map_or(false, |ext| ext == "gz"))
        }

        pub fn writer(&self) -> Result<Box<dyn Write>, Error> {
            match &self.filepath {
                Some(path) => {
                    let file = File::create(path)?;
                    let writer: Box<dyn Write> = if self.should_compress() {
                        Box::new(BufWriter::with_capacity(
                            self.buffer_size,
                            GzEncoder::new(file, self.compression_level),
                        ))
                    } else {
                        Box::new(BufWriter::with_capacity(self.buffer_size, file))
                    };
                    Ok(writer)
                }
                None => Ok(Box::new(BufWriter::with_capacity(
                    self.buffer_size,
                    io::stdout(),
                ))),
            }
        }
    }
}
