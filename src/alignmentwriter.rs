use std::{ffi::CString, io, ptr::NonNull};

use rust_htslib::{bam, htslib};

pub(crate) struct AlignmentWriter {
    file: Option<NonNull<htslib::htsFile>>,
    header: bam::HeaderView,
}

impl AlignmentWriter {
    pub fn new(path: &str, header: &bam::Header, binary: bool) -> io::Result<Self> {
        let path = CString::new(path)
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidInput, error))?;
        let mode: &[u8] = if binary { b"wb\0" } else { b"w\0" };
        let header = bam::HeaderView::from_header(header);
        // Both strings are NUL-terminated and live for the duration of the call.
        let file = unsafe { htslib::hts_open(path.as_ptr(), mode.as_ptr().cast()) };
        let writer = Self {
            file: Some(NonNull::new(file).ok_or_else(io::Error::last_os_error)?),
            header,
        };
        // The file is open and this writer owns the header for its entire lifetime.
        if unsafe { htslib::sam_hdr_write(writer.file.unwrap().as_ptr(), writer.header.inner()) }
            < 0
        {
            return Err(io::Error::other("Failed to write alignment header"));
        }
        Ok(writer)
    }

    pub fn write(&mut self, record: &bam::Record) -> io::Result<()> {
        // Header and record remain alive while HTSlib serializes this record.
        if unsafe {
            htslib::sam_write1(
                self.file.unwrap().as_ptr(),
                self.header.inner(),
                record.inner(),
            )
        } < 0
        {
            return Err(io::Error::other("Failed to write alignment record"));
        }
        Ok(())
    }

    pub fn finish(mut self) -> io::Result<()> {
        let file = self.file.take().unwrap();
        // Taking the sole handle prevents Drop from closing it a second time.
        if unsafe { htslib::hts_close(file.as_ptr()) } != 0 {
            return Err(io::Error::other("Failed to finalize alignment output"));
        }
        Ok(())
    }
}

impl Drop for AlignmentWriter {
    fn drop(&mut self) {
        if let Some(file) = self.file.take() {
            // Best-effort cleanup during errors; successful runs must call finish.
            unsafe {
                htslib::hts_close(file.as_ptr());
            }
        }
    }
}
