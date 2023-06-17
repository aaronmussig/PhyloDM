#[derive(Debug)]
pub struct PhyloErr(pub String);

impl std::fmt::Display for PhyloErr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for PhyloErr {}
