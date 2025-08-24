use std::path::Path;

pub struct BenchmarkData {
    pub bam_path: String,
}

impl BenchmarkData {
    pub fn new(bam_path: String) -> Self {
        Self { bam_path }
    }
}

pub fn generate_small_dataset() -> Result<BenchmarkData, Box<dyn std::error::Error>> {
    if Path::new("tests/test1.bam").exists() {
        Ok(BenchmarkData::new("tests/test1.bam".to_string()))
    } else {
        Err("Test BAM file not found. Run tests first.".into())
    }
}

pub fn generate_medium_dataset() -> Result<BenchmarkData, Box<dyn std::error::Error>> {
    if Path::new("tests/test2.bam").exists() {
        Ok(BenchmarkData::new("tests/test2.bam".to_string()))
    } else if Path::new("tests/test1.bam").exists() {
        Ok(BenchmarkData::new("tests/test1.bam".to_string()))
    } else {
        Err("Test BAM file not found. Run tests first.".into())
    }
}

pub fn generate_large_dataset() -> Result<BenchmarkData, Box<dyn std::error::Error>> {
    if Path::new("tests/test3.bam").exists() {
        Ok(BenchmarkData::new("tests/test3.bam".to_string()))
    } else if Path::new("tests/test1.bam").exists() {
        Ok(BenchmarkData::new("tests/test1.bam".to_string()))
    } else {
        Err("Test BAM file not found. Run tests first.".into())
    }
}

pub fn generate_high_depth_dataset() -> Result<BenchmarkData, Box<dyn std::error::Error>> {
    if Path::new("tests/test4.bam").exists() {
        Ok(BenchmarkData::new("tests/test4.bam".to_string()))
    } else if Path::new("tests/test1.bam").exists() {
        Ok(BenchmarkData::new("tests/test1.bam".to_string()))
    } else {
        Err("Test BAM file not found. Run tests first.".into())
    }
}

pub fn generate_low_methylation_dataset() -> Result<BenchmarkData, Box<dyn std::error::Error>> {
    if Path::new("tests/test5.bam").exists() {
        Ok(BenchmarkData::new("tests/test5.bam".to_string()))
    } else if Path::new("tests/test1.bam").exists() {
        Ok(BenchmarkData::new("tests/test1.bam".to_string()))
    } else {
        Err("Test BAM file not found. Run tests first.".into())
    }
}

pub fn generate_high_discordance_dataset() -> Result<BenchmarkData, Box<dyn std::error::Error>> {
    if Path::new("tests/test6.bam").exists() {
        Ok(BenchmarkData::new("tests/test6.bam".to_string()))
    } else if Path::new("tests/test1.bam").exists() {
        Ok(BenchmarkData::new("tests/test1.bam".to_string()))
    } else {
        Err("Test BAM file not found. Run tests first.".into())
    }
}