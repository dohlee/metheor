use std::collections::{BTreeSet, HashMap, HashSet};

use rust_htslib::bam::{record::Aux, Record};

use crate::readutil::{BismarkRead, CpGPosition, Quartet};

type ReadKey = (Vec<u8>, Option<Vec<u8>>);

#[derive(Default)]
pub(crate) struct Observations {
    pub cpgs: HashSet<CpGPosition>,
    #[allow(dead_code)] // Quartet metrics are compiled only in the binary crate.
    pub quartets: HashSet<Quartet>,
}

struct PendingMate {
    position: CpGPosition,
    mate: CpGPosition,
    first_in_template: bool,
    observations: Observations,
}

pub(crate) struct OverlapTracker {
    enabled: bool,
    pending: HashMap<ReadKey, PendingMate>,
    deadlines: BTreeSet<(CpGPosition, ReadKey)>,
    last_position: Option<CpGPosition>,
}

impl OverlapTracker {
    pub fn new(enabled: bool) -> Self {
        Self {
            enabled,
            pending: HashMap::new(),
            deadlines: BTreeSet::new(),
            last_position: None,
        }
    }

    pub fn observe(
        &mut self,
        record: &Record,
        read: &BismarkRead,
        counted_cpgs: &[CpGPosition],
    ) -> Observations {
        if !self.enabled || record.is_unmapped() {
            return Observations::default();
        }

        let position = CpGPosition::new(record.tid(), record.pos() as i32);
        if let Some(previous) = self.last_position {
            assert!(
                previous <= position,
                "--ignore-overlaps requires coordinate-sorted alignments"
            );
        }
        self.last_position = Some(position);

        // A missing mate must not keep its observation set for the rest of the file.
        while self
            .deadlines
            .first()
            .is_some_and(|(deadline, _)| *deadline < position)
        {
            let (deadline, key) = self.deadlines.pop_first().unwrap();
            if self
                .pending
                .get(&key)
                .is_some_and(|pending| pending.mate == deadline)
            {
                self.pending.remove(&key);
            }
        }

        if !record.is_paired()
            || record.qname() == b"*"
            || record.is_mate_unmapped()
            || record.is_secondary()
            || record.is_supplementary()
            || record.is_first_in_template() == record.is_last_in_template()
            || record.tid() != record.mtid()
            || record.mpos() < 0
        {
            return Observations::default();
        }

        let group = match record.aux(b"RG") {
            Ok(Aux::String(group)) => Some(group.as_bytes().to_vec()),
            _ => None,
        };
        let key = (record.qname().to_vec(), group);
        let mate = CpGPosition::new(record.mtid(), record.mpos() as i32);

        if self.pending.get(&key).is_some_and(|pending| {
            pending.position == mate
                && pending.mate == position
                && pending.first_in_template != record.is_first_in_template()
        }) {
            let pending = self.pending.remove(&key).unwrap();
            self.deadlines.remove(&(position, key));
            return pending.observations;
        }

        if mate < position || mate.pos > read.get_end_pos() + 1 {
            return Observations::default();
        }

        self.deadlines.insert((mate, key.clone()));
        self.pending.insert(
            key,
            PendingMate {
                position,
                mate,
                first_in_template: record.is_first_in_template(),
                observations: Observations {
                    cpgs: counted_cpgs.iter().copied().collect(),
                    quartets: read.get_cpg_quartets_and_patterns().0.into_iter().collect(),
                },
            },
        );
        Observations::default()
    }
}
