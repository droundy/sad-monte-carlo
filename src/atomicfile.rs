use tempfile::tempfile_in;

use std::path;
use std::fs::File;
use std::io::{Result};

pub struct AtomicFile {
    path: path::PathBuf,
    file: File,
}

impl AtomicFile {
    fn create<P: AsPath>(p: P) -> Result<AtomicFile> {
        Ok(AtomicFile {
            path: path::PathBuf::from(p),
            file: File::create(p),
        })
    }
}

fn safe_parent(p: &path::Path) -> Option<&path::Path> {
    match p.parent() {
        None => None,
        Some(x) if x.as_os_str().len() == 0 => Some(&path::Path::new(".")),
        x => x
    }
}

