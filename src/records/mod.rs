// records/mod.rs
use std::fmt;

pub trait Record: Sized + for<'a> From<Self::Slice<'a>> {
    type Slice<'a>: RecordSlice<'a, Owned = Self>;
    fn start(&self) -> u32;
    fn end(&self) -> u32;
    fn to_bytes(&self) -> Vec<u8>;
}

pub trait RecordSlice<'a>: Sized {
    type Owned: Record + From<Self>;
    fn from_bytes(bytes: &'a [u8]) -> Self;
    fn start(&self) -> u32;
    fn end(&self) -> u32;
    fn to_owned(self) -> Self::Owned;
}

#[derive(Debug, Clone, PartialEq)]
pub struct BedRecord {
    pub start: u32,
    pub end: u32,
    pub rest: String,
}

#[derive(Debug, PartialEq)]
pub struct BedRecordSlice<'a> {
    pub start: u32,
    pub end: u32,
    pub rest: &'a [u8],
}

impl Record for BedRecord {
    type Slice<'a> = BedRecordSlice<'a>;

    fn start(&self) -> u32 {
        self.start
    }
    fn end(&self) -> u32 {
        self.end
    }

    fn to_bytes(&self) -> Vec<u8> {
        // manual serialization
        let mut bytes = Vec::with_capacity(8 + self.rest.len());
        bytes.extend_from_slice(&self.start.to_le_bytes());
        bytes.extend_from_slice(&self.end.to_le_bytes());
        bytes.extend_from_slice(self.rest.as_bytes());
        bytes
    }
}

impl<'a> RecordSlice<'a> for BedRecordSlice<'a> {
    type Owned = BedRecord;

    fn start(&self) -> u32 {
        self.start
    }
    fn end(&self) -> u32 {
        self.end
    }

    fn from_bytes(bytes: &'a [u8]) -> Self {
        if bytes.len() < 8 {
            panic!("Internal error: invalid byte record, bytes length too small.")
        }

        // SAFETY: We've checked the length above, and we know the slices are 4 bytes each
        // u32.
        //unsafe {
        //    Self {
        //        start: u32::from_le_bytes(*(bytes.as_ptr() as *const [u8; 4])),
        //        end: u32::from_le_bytes(*(bytes[4..].as_ptr() as *const [u8; 4])),
        //        rest: &bytes[8..],
        //    }
        //}
        unsafe {
            let start = u32::from_le_bytes(*(bytes.get_unchecked(0..4).as_ptr() as *const [u8; 4]));
            let end = u32::from_le_bytes(*(bytes.get_unchecked(4..8).as_ptr() as *const [u8; 4]));
            let rest = bytes.get_unchecked(8..);
            Self { start, end, rest }
        }
    }

    fn to_owned(self) -> Self::Owned {
        BedRecord {
            start: self.start,
            end: self.end,
            rest: std::str::from_utf8(self.rest).unwrap().to_string(),
        }
    }
}

impl From<BedRecordSlice<'_>> for BedRecord {
    fn from(slice: BedRecordSlice<'_>) -> Self {
        Self {
            start: slice.start,
            end: slice.end,
            rest: std::str::from_utf8(slice.rest).unwrap().to_string(),
        }
    }
}

impl fmt::Display for BedRecordSlice<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.rest.is_empty() {
            write!(f, "{}\t{}", self.start, self.end)
        } else {
            // Assume ASCII and use lossy UTF8 conversion for display only
            write!(
                f,
                "{}\t{}\t{}",
                self.start,
                self.end,
                String::from_utf8_lossy(self.rest)
            )
        }
    }
}

// // Just use derive(Debug) instead of manual impls
// impl fmt::Display for BedRecord {
//     fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
//         if self.rest.is_empty() {
//             write!(f, "{}\t{}", self.start, self.end)
//         } else {
//             write!(f, "{}\t{}\t{}", self.start, self.end, self.rest)
//         }
//     }
// }
//
// impl fmt::Display for BedRecordSlice<'_> {
//     fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
//         if self.rest.is_empty() {
//             write!(f, "{}\t{}", self.start, self.end)
//         } else {
//             write!(f, "{}\t{}\t{}", self.start, self.end, self.rest)
//         }
//     }
// }
