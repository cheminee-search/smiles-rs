//! Opt-in, versioned preparation contracts for RDKit 2026.09.1.
//!
//! These APIs do not change ordinary [`crate::Smiles`] parsing or its legacy
//! radical policy. Each capability documents its own admitted graph domain;
//! this module does not provide a general native sanitization implementation.

mod element_data;
mod radicals;

pub use radicals::{
    AtomRadicalState, NativeValenceCache, RadicalDecision, RadicalPreparationError,
    RadicalPreparedSmiles, UnsupportedPreparation,
};
