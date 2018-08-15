//! A simple thing to atomically write to a file.

use tempfile::TempDir;

use std::path::{PathBuf, Path};
use std::fs::{File, rename};
use std::io::{Result, Error, ErrorKind, Write};

/// A version of File that should never leave a partially-written
/// file.  This is only useful for creating files, and will overwrite
/// an existing file with the same name.
pub struct AtomicFile {
    path: PathBuf,
    dir: TempDir,
    file: File,
}

impl AtomicFile {
    /// Create a file.
    pub fn create<P: AsRef<Path>>(p: P) -> Result<AtomicFile> {
        let filepath = p.as_ref();
        let p = match filepath.parent() {
            None => {
                return Err(Error::new(ErrorKind::Other,
                                      format!("Cannot create a file named {:?}",
                                              filepath)));
            }
            Some(p) if p.as_os_str().len() == 0 => {
                &Path::new(".")
            }
            Some(p) => p,
        };
        let dir = TempDir::new_in(p)?;
        let file_path = dir.path().join("temp");
        let file = File::create(file_path)?;
        Ok(AtomicFile {
            path: PathBuf::from(filepath),
            dir: dir,
            file: file,
        })
    }
}

impl<'a> Write for &'a AtomicFile {
    fn write(&mut self, buf: &[u8]) -> Result<usize> {
        (&self.file).write(buf)
    }
    fn flush(&mut self) -> Result<()> {
        (&self.file).flush()
    }

    fn write_all(&mut self, buf: &[u8]) -> Result<()> {
        (&self.file).write_all(buf)
    }
}

impl Write for AtomicFile {
    fn write(&mut self, buf: &[u8]) -> Result<usize> {
        self.file.write(buf)
    }
    fn flush(&mut self) -> Result<()> {
        self.file.flush()
    }

    fn write_all(&mut self, buf: &[u8]) -> Result<()> {
        self.file.write_all(buf)
    }
}

impl Drop for AtomicFile {
    fn drop(&mut self) {
        self.file.sync_data().ok();
        // we ignore errors in drop, because there is no nice way to
        // handle them.
        rename( self.dir.path().join("temp"), &self.path ).ok();
    }
}
