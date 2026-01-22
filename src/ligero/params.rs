//! Ligero parameters.
//!
//! This module defines the parameters for a Ligero commitment scheme instance.

/// Parameters for a Ligero commitment scheme.
#[derive(Clone, Debug)]
pub struct LigeroParams {
    /// Number of columns requested by verifier (random padding width).
    pub nreq: usize,

    /// Inverse rate of the error correcting code.
    pub rate: usize,

    /// Size of each row in terms of field elements (BLOCK).
    pub block: usize,

    /// Extended block size (DBLOCK = 2 * BLOCK - 1).
    pub dblock: usize,

    /// Number of witness values per row.
    pub wr: usize,

    /// Number of quadratic constraints per row triple.
    pub qr: usize,

    /// Row index where witness values start (typically 3).
    pub iw: usize,

    /// Row index where quadratic constraints begin.
    pub iq: usize,

    /// Number of linear constraints.
    pub nl: usize,

    /// Number of quadratic constraints.
    pub nq: usize,

    /// Number of rows used for witnesses.
    pub nwrow: usize,

    /// Number of row triples for quadratic constraints.
    pub nqt: usize,

    /// Total rows for witnesses and quadratic constraints (NWROW + NQT).
    pub nqw: usize,

    /// Total number of rows in the tableau.
    pub nrow: usize,

    /// Total number of columns in the tableau.
    pub ncol: usize,
}

impl LigeroParams {
    /// Create Ligero parameters for a given witness size and number of quadratic constraints.
    ///
    /// # Arguments
    /// * `witness_size` - The number of witness elements
    /// * `num_quadratic` - The number of quadratic constraints
    /// * `security_level` - The desired security level (typically 128)
    pub fn new(witness_size: usize, num_quadratic: usize, security_level: usize) -> Self {
        // Determine NREQ based on security level
        // For 128-bit security with rate 4, we need about 6-10 column queries
        let nreq = match security_level {
            128 => 6,
            256 => 12,
            _ => ((security_level + 20) / 21).max(4),
        };

        // Rate is typically 4 for good efficiency/soundness tradeoff
        let rate = 4;

        // WR and QR must satisfy WR >= QR >= NREQ
        let qr = nreq.max(2);
        let wr = qr.max(witness_size.min(20));

        // BLOCK = NREQ + WR
        let block = nreq + wr;

        // Ensure BLOCK >= 2 * (NREQ + WR) - 1 for degree-2 polynomial product
        let block = block.max(2 * (nreq + wr) - 1);

        let dblock = 2 * block - 1;

        // Number of witness rows
        let nwrow = (witness_size + wr - 1) / wr;

        // Number of quadratic constraint triples
        let nqt = (num_quadratic + qr - 1) / qr;

        // Starting indices
        let iw = 3; // After ILDT, IDOT, IQD
        let iq = iw + nwrow;

        // Total rows: 3 ZK rows + witness rows + 3 * quadratic triples
        let nrow = 3 + nwrow + 3 * nqt;

        // NQW: witness rows + quadratic rows (one per triple, not 3)
        let nqw = nwrow + nqt;

        // Total columns
        let ncol = block * rate;

        Self {
            nreq,
            rate,
            block,
            dblock,
            wr,
            qr,
            iw,
            iq,
            nl: 0, // Set by the prover based on linear constraints
            nq: num_quadratic,
            nwrow,
            nqt,
            nqw,
            nrow,
            ncol,
        }
    }

    /// Create parameters from explicit values (for test vectors).
    pub fn from_explicit(
        nreq: usize,
        rate: usize,
        wr: usize,
        qr: usize,
        nrow: usize,
        nq: usize,
    ) -> Self {
        let block = nreq + wr;
        let dblock = 2 * block - 1;
        let ncol = block * rate;

        // Compute derived values
        let nqt = (nq + qr - 1) / qr.max(1);
        let nwrow = nrow - 3 - 3 * nqt;
        let iw = 3;
        let iq = iw + nwrow;
        let nqw = nwrow + nqt;

        Self {
            nreq,
            rate,
            block,
            dblock,
            wr,
            qr,
            iw,
            iq,
            nl: 0,
            nq,
            nwrow,
            nqt,
            nqw,
            nrow,
            ncol,
        }
    }

    /// Get the row index for the low-degree test row.
    pub const fn ildt(&self) -> usize {
        0
    }

    /// Get the row index for the dot product test row.
    pub const fn idot(&self) -> usize {
        1
    }

    /// Get the row index for the quadratic test row.
    pub const fn iqd(&self) -> usize {
        2
    }

    /// Get the row index for the i-th witness row.
    pub fn witness_row(&self, i: usize) -> usize {
        self.iw + i
    }

    /// Get the row index for the x-component of the i-th quadratic triple.
    pub fn quadratic_x_row(&self, i: usize) -> usize {
        self.iq + 3 * i
    }

    /// Get the row index for the y-component of the i-th quadratic triple.
    pub fn quadratic_y_row(&self, i: usize) -> usize {
        self.iq + 3 * i + 1
    }

    /// Get the row index for the z-component of the i-th quadratic triple.
    pub fn quadratic_z_row(&self, i: usize) -> usize {
        self.iq + 3 * i + 2
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_params_creation() {
        let params = LigeroParams::new(100, 10, 128);

        // Basic sanity checks
        assert!(params.block >= params.nreq + params.wr);
        assert!(params.ncol >= params.block);
        assert!(params.nrow >= 3);
        assert!(params.wr >= params.qr);
    }

    #[test]
    fn test_explicit_params() {
        // From test vector: NREQ=6, RATE=4, WR=20, QR=2, NROW=7, NQ=1
        let params = LigeroParams::from_explicit(6, 4, 20, 2, 7, 1);

        assert_eq!(params.nreq, 6);
        assert_eq!(params.rate, 4);
        assert_eq!(params.wr, 20);
        assert_eq!(params.qr, 2);
        assert_eq!(params.nrow, 7);
        assert_eq!(params.nq, 1);
        assert_eq!(params.block, 26); // NREQ + WR = 6 + 20
    }
}
