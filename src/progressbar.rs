pub struct ProgressBar {
    bar: indicatif::ProgressBar,
}

impl Default for ProgressBar {
    fn default() -> Self {
        Self::new()
    }
}

impl ProgressBar {
    pub fn new() -> Self {
        let bar = indicatif::ProgressBar::new(1);
        bar.set_style(
            indicatif::ProgressStyle::default_bar()
                .template("{spinner} {elapsed_precise} {msg}")
                .expect("Invalid progress bar template"),
        );

        Self { bar }
    }

    pub fn inc_length(&self, i: u64) {
        self.bar.inc_length(i);
    }

    pub fn inc(&self, i: u64) {
        self.bar.inc(i);
    }

    pub fn set_message(&self, s: String) {
        self.bar.set_message(s);
    }

    pub fn update(&self, readcount: i32, valid_readcount: i32) {
        self.inc_length(10000);
        self.inc(10000);
        self.set_message(format!(
            "Processed {} reads, found {} valid reads.",
            readcount, valid_readcount
        ));
    }

    pub fn update_lpmd(&self, progress_string: String) {
        self.inc_length(10000);
        self.inc(10000);
        self.set_message(progress_string);
    }
}
