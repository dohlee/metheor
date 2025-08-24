use std::time::{Duration, Instant};

#[derive(Debug, Clone)]
pub struct PerformanceMetrics {
    pub measure: String,
    pub dataset_size: usize,
    pub execution_time: Duration,
    pub reads_per_second: f64,
    pub memory_usage_mb: Option<f64>,
}

pub struct PerformanceAnalyzer {
    metrics: Vec<PerformanceMetrics>,
}

impl PerformanceAnalyzer {
    pub fn new() -> Self {
        Self {
            metrics: Vec::new(),
        }
    }

    pub fn record_metric(&mut self, metric: PerformanceMetrics) {
        self.metrics.push(metric);
    }

    pub fn analyze_scaling(&self, measure: &str) -> ScalingAnalysis {
        let mut measure_metrics: Vec<_> = self.metrics
            .iter()
            .filter(|m| m.measure == measure)
            .cloned()
            .collect();
        
        measure_metrics.sort_by_key(|m| m.dataset_size);
        
        let scaling_factor = if measure_metrics.len() >= 2 {
            let first = &measure_metrics[0];
            let last = &measure_metrics[measure_metrics.len() - 1];
            
            let size_ratio = last.dataset_size as f64 / first.dataset_size as f64;
            let time_ratio = last.execution_time.as_secs_f64() / first.execution_time.as_secs_f64();
            
            time_ratio / size_ratio
        } else {
            1.0
        };
        
        ScalingAnalysis {
            measure: measure.to_string(),
            scaling_factor,
            is_linear: (0.9..=1.1).contains(&scaling_factor),
            metrics: measure_metrics,
        }
    }

    pub fn compare_measures(&self) -> Vec<MeasureComparison> {
        let mut comparisons = Vec::new();
        let measures: Vec<String> = self.metrics
            .iter()
            .map(|m| m.measure.clone())
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();
        
        for measure in measures {
            let avg_time = self.average_execution_time(&measure);
            let avg_throughput = self.average_throughput(&measure);
            
            comparisons.push(MeasureComparison {
                measure,
                average_time: avg_time,
                average_throughput: avg_throughput,
            });
        }
        
        comparisons.sort_by(|a, b| a.average_time.partial_cmp(&b.average_time).unwrap());
        comparisons
    }

    fn average_execution_time(&self, measure: &str) -> Duration {
        let times: Vec<Duration> = self.metrics
            .iter()
            .filter(|m| m.measure == measure)
            .map(|m| m.execution_time)
            .collect();
        
        if times.is_empty() {
            return Duration::from_secs(0);
        }
        
        let total: Duration = times.iter().sum();
        total / times.len() as u32
    }

    fn average_throughput(&self, measure: &str) -> f64 {
        let throughputs: Vec<f64> = self.metrics
            .iter()
            .filter(|m| m.measure == measure)
            .map(|m| m.reads_per_second)
            .collect();
        
        if throughputs.is_empty() {
            return 0.0;
        }
        
        throughputs.iter().sum::<f64>() / throughputs.len() as f64
    }

    pub fn generate_report(&self) -> String {
        let mut report = String::new();
        
        report.push_str("# Metheor Performance Analysis Report\n\n");
        
        report.push_str("## Performance Summary\n\n");
        let comparisons = self.compare_measures();
        report.push_str("| Measure | Avg Time (ms) | Avg Throughput (reads/s) |\n");
        report.push_str("|---------|---------------|-------------------------|\n");
        for comp in comparisons {
            report.push_str(&format!(
                "| {} | {:.2} | {:.0} |\n",
                comp.measure,
                comp.average_time.as_millis(),
                comp.average_throughput
            ));
        }
        
        report.push_str("\n## Scaling Analysis\n\n");
        let measures: Vec<String> = self.metrics
            .iter()
            .map(|m| m.measure.clone())
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();
        
        for measure in measures {
            let analysis = self.analyze_scaling(&measure);
            report.push_str(&format!("### {}\n", measure));
            report.push_str(&format!(
                "- Scaling Factor: {:.2}x\n",
                analysis.scaling_factor
            ));
            report.push_str(&format!(
                "- Scaling Type: {}\n",
                if analysis.is_linear { "Linear" } else { "Non-linear" }
            ));
            report.push_str("\n");
        }
        
        report
    }
}

#[derive(Debug)]
pub struct ScalingAnalysis {
    pub measure: String,
    pub scaling_factor: f64,
    pub is_linear: bool,
    pub metrics: Vec<PerformanceMetrics>,
}

#[derive(Debug)]
pub struct MeasureComparison {
    pub measure: String,
    pub average_time: Duration,
    pub average_throughput: f64,
}

pub struct BenchmarkProfiler {
    start_time: Option<Instant>,
    measure_name: String,
}

impl BenchmarkProfiler {
    pub fn new(measure_name: &str) -> Self {
        Self {
            start_time: None,
            measure_name: measure_name.to_string(),
        }
    }

    pub fn start(&mut self) {
        self.start_time = Some(Instant::now());
    }

    pub fn stop(&self, dataset_size: usize) -> PerformanceMetrics {
        let duration = self.start_time
            .map(|start| start.elapsed())
            .unwrap_or_else(|| Duration::from_secs(0));
        
        let reads_per_second = if duration.as_secs_f64() > 0.0 {
            dataset_size as f64 / duration.as_secs_f64()
        } else {
            0.0
        };
        
        PerformanceMetrics {
            measure: self.measure_name.clone(),
            dataset_size,
            execution_time: duration,
            reads_per_second,
            memory_usage_mb: None,
        }
    }
}

pub fn estimate_memory_usage() -> Option<f64> {
    #[cfg(target_os = "linux")]
    {
        if let Ok(status) = std::fs::read_to_string("/proc/self/status") {
            for line in status.lines() {
                if line.starts_with("VmRSS:") {
                    let parts: Vec<&str> = line.split_whitespace().collect();
                    if parts.len() >= 2 {
                        if let Ok(kb) = parts[1].parse::<f64>() {
                            return Some(kb / 1024.0);
                        }
                    }
                }
            }
        }
    }
    
    None
}

pub fn format_duration(duration: Duration) -> String {
    let millis = duration.as_millis();
    if millis < 1000 {
        format!("{}ms", millis)
    } else {
        format!("{:.2}s", duration.as_secs_f64())
    }
}

pub fn format_throughput(reads_per_second: f64) -> String {
    if reads_per_second > 1_000_000.0 {
        format!("{:.2}M reads/s", reads_per_second / 1_000_000.0)
    } else if reads_per_second > 1_000.0 {
        format!("{:.2}K reads/s", reads_per_second / 1_000.0)
    } else {
        format!("{:.0} reads/s", reads_per_second)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_performance_analyzer() {
        let mut analyzer = PerformanceAnalyzer::new();
        
        analyzer.record_metric(PerformanceMetrics {
            measure: "PDR".to_string(),
            dataset_size: 1000,
            execution_time: Duration::from_millis(100),
            reads_per_second: 10000.0,
            memory_usage_mb: Some(50.0),
        });
        
        analyzer.record_metric(PerformanceMetrics {
            measure: "PDR".to_string(),
            dataset_size: 10000,
            execution_time: Duration::from_millis(1000),
            reads_per_second: 10000.0,
            memory_usage_mb: Some(100.0),
        });
        
        let scaling = analyzer.analyze_scaling("PDR");
        assert!(scaling.is_linear);
    }
    
    #[test]
    fn test_format_duration() {
        assert_eq!(format_duration(Duration::from_millis(500)), "500ms");
        assert_eq!(format_duration(Duration::from_millis(1500)), "1.50s");
    }
    
    #[test]
    fn test_format_throughput() {
        assert_eq!(format_throughput(500.0), "500 reads/s");
        assert_eq!(format_throughput(5_000.0), "5.00K reads/s");
        assert_eq!(format_throughput(5_000_000.0), "5.00M reads/s");
    }
}