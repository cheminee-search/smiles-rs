//! Native element-table data for the versioned admission rule.
//!
//! Derived from RDKit atomic_data.cpp at commit
//! fece8caa860bdf6c9c82bdb9253f22a0cf87150c. The source BSD-3-Clause
//! notice is retained in LICENSES/RDKit-BSD-3-Clause.txt.

#[derive(Clone, Copy)]
pub(super) struct ElementData {
    pub(super) outer: u8,
    pub(super) valences: &'static [i8],
}

pub(super) const ELEMENTS: [ElementData; 119] = [
    ElementData { outer: 0, valences: &[-1] },         // 0 *
    ElementData { outer: 1, valences: &[1] },          // 1 H
    ElementData { outer: 2, valences: &[0] },          // 2 He
    ElementData { outer: 1, valences: &[1, -1] },      // 3 Li
    ElementData { outer: 2, valences: &[2] },          // 4 Be
    ElementData { outer: 3, valences: &[3] },          // 5 B
    ElementData { outer: 4, valences: &[4] },          // 6 C
    ElementData { outer: 5, valences: &[3] },          // 7 N
    ElementData { outer: 6, valences: &[2] },          // 8 O
    ElementData { outer: 7, valences: &[1] },          // 9 F
    ElementData { outer: 8, valences: &[0] },          // 10 Ne
    ElementData { outer: 1, valences: &[1, -1] },      // 11 Na
    ElementData { outer: 2, valences: &[2, -1] },      // 12 Mg
    ElementData { outer: 3, valences: &[3] },          // 13 Al
    ElementData { outer: 4, valences: &[4] },          // 14 Si
    ElementData { outer: 5, valences: &[3, 5] },       // 15 P
    ElementData { outer: 6, valences: &[2, 4, 6] },    // 16 S
    ElementData { outer: 7, valences: &[1] },          // 17 Cl
    ElementData { outer: 8, valences: &[0] },          // 18 Ar
    ElementData { outer: 1, valences: &[1, -1] },      // 19 K
    ElementData { outer: 2, valences: &[2, -1] },      // 20 Ca
    ElementData { outer: 3, valences: &[-1] },         // 21 Sc
    ElementData { outer: 4, valences: &[-1] },         // 22 Ti
    ElementData { outer: 5, valences: &[-1] },         // 23 V
    ElementData { outer: 6, valences: &[-1] },         // 24 Cr
    ElementData { outer: 7, valences: &[-1] },         // 25 Mn
    ElementData { outer: 8, valences: &[-1] },         // 26 Fe
    ElementData { outer: 9, valences: &[-1] },         // 27 Co
    ElementData { outer: 10, valences: &[-1] },        // 28 Ni
    ElementData { outer: 11, valences: &[-1] },        // 29 Cu
    ElementData { outer: 2, valences: &[-1] },         // 30 Zn
    ElementData { outer: 3, valences: &[3] },          // 31 Ga
    ElementData { outer: 4, valences: &[4] },          // 32 Ge
    ElementData { outer: 5, valences: &[3, 5] },       // 33 As
    ElementData { outer: 6, valences: &[2, 4, 6] },    // 34 Se
    ElementData { outer: 7, valences: &[1] },          // 35 Br
    ElementData { outer: 8, valences: &[0] },          // 36 Kr
    ElementData { outer: 1, valences: &[1, -1] },      // 37 Rb
    ElementData { outer: 2, valences: &[2, -1] },      // 38 Sr
    ElementData { outer: 3, valences: &[-1] },         // 39 Y
    ElementData { outer: 4, valences: &[-1] },         // 40 Zr
    ElementData { outer: 5, valences: &[-1] },         // 41 Nb
    ElementData { outer: 6, valences: &[-1] },         // 42 Mo
    ElementData { outer: 7, valences: &[-1] },         // 43 Tc
    ElementData { outer: 8, valences: &[-1] },         // 44 Ru
    ElementData { outer: 9, valences: &[-1] },         // 45 Rh
    ElementData { outer: 10, valences: &[-1] },        // 46 Pd
    ElementData { outer: 11, valences: &[-1] },        // 47 Ag
    ElementData { outer: 2, valences: &[-1] },         // 48 Cd
    ElementData { outer: 3, valences: &[3] },          // 49 In
    ElementData { outer: 4, valences: &[2, 4] },       // 50 Sn
    ElementData { outer: 5, valences: &[3, 5] },       // 51 Sb
    ElementData { outer: 6, valences: &[2, 4, 6] },    // 52 Te
    ElementData { outer: 7, valences: &[1, 3, 5] },    // 53 I
    ElementData { outer: 8, valences: &[0, 2, 4, 6] }, // 54 Xe
    ElementData { outer: 1, valences: &[1] },          // 55 Cs
    ElementData { outer: 2, valences: &[2, -1] },      // 56 Ba
    ElementData { outer: 3, valences: &[-1] },         // 57 La
    ElementData { outer: 4, valences: &[-1] },         // 58 Ce
    ElementData { outer: 3, valences: &[-1] },         // 59 Pr
    ElementData { outer: 4, valences: &[-1] },         // 60 Nd
    ElementData { outer: 5, valences: &[-1] },         // 61 Pm
    ElementData { outer: 6, valences: &[-1] },         // 62 Sm
    ElementData { outer: 7, valences: &[-1] },         // 63 Eu
    ElementData { outer: 8, valences: &[-1] },         // 64 Gd
    ElementData { outer: 9, valences: &[-1] },         // 65 Tb
    ElementData { outer: 10, valences: &[-1] },        // 66 Dy
    ElementData { outer: 11, valences: &[-1] },        // 67 Ho
    ElementData { outer: 12, valences: &[-1] },        // 68 Er
    ElementData { outer: 13, valences: &[-1] },        // 69 Tm
    ElementData { outer: 14, valences: &[-1] },        // 70 Yb
    ElementData { outer: 15, valences: &[-1] },        // 71 Lu
    ElementData { outer: 4, valences: &[-1] },         // 72 Hf
    ElementData { outer: 5, valences: &[-1] },         // 73 Ta
    ElementData { outer: 6, valences: &[-1] },         // 74 W
    ElementData { outer: 7, valences: &[-1] },         // 75 Re
    ElementData { outer: 8, valences: &[-1] },         // 76 Os
    ElementData { outer: 9, valences: &[-1] },         // 77 Ir
    ElementData { outer: 10, valences: &[-1] },        // 78 Pt
    ElementData { outer: 11, valences: &[-1] },        // 79 Au
    ElementData { outer: 2, valences: &[-1] },         // 80 Hg
    ElementData { outer: 3, valences: &[-1] },         // 81 Tl
    ElementData { outer: 4, valences: &[2, 4] },       // 82 Pb
    ElementData { outer: 5, valences: &[3, 5] },       // 83 Bi
    ElementData { outer: 6, valences: &[2, 4, 6] },    // 84 Po
    ElementData { outer: 7, valences: &[1, 3, 5] },    // 85 At
    ElementData { outer: 8, valences: &[0] },          // 86 Rn
    ElementData { outer: 1, valences: &[1] },          // 87 Fr
    ElementData { outer: 2, valences: &[2, -1] },      // 88 Ra
    ElementData { outer: 3, valences: &[-1] },         // 89 Ac
    ElementData { outer: 4, valences: &[-1] },         // 90 Th
    ElementData { outer: 3, valences: &[-1] },         // 91 Pa
    ElementData { outer: 4, valences: &[-1] },         // 92 U
    ElementData { outer: 5, valences: &[-1] },         // 93 Np
    ElementData { outer: 6, valences: &[-1] },         // 94 Pu
    ElementData { outer: 7, valences: &[-1] },         // 95 Am
    ElementData { outer: 8, valences: &[-1] },         // 96 Cm
    ElementData { outer: 9, valences: &[-1] },         // 97 Bk
    ElementData { outer: 10, valences: &[-1] },        // 98 Cf
    ElementData { outer: 11, valences: &[-1] },        // 99 Es
    ElementData { outer: 12, valences: &[-1] },        // 100 Fm
    ElementData { outer: 13, valences: &[-1] },        // 101 Md
    ElementData { outer: 14, valences: &[-1] },        // 102 No
    ElementData { outer: 15, valences: &[-1] },        // 103 Lr
    ElementData { outer: 2, valences: &[-1] },         // 104 Rf
    ElementData { outer: 2, valences: &[-1] },         // 105 Db
    ElementData { outer: 2, valences: &[-1] },         // 106 Sg
    ElementData { outer: 2, valences: &[-1] },         // 107 Bh
    ElementData { outer: 2, valences: &[-1] },         // 108 Hs
    ElementData { outer: 2, valences: &[-1] },         // 109 Mt
    ElementData { outer: 2, valences: &[-1] },         // 110 Ds
    ElementData { outer: 2, valences: &[-1] },         // 111 Rg
    ElementData { outer: 2, valences: &[-1] },         // 112 Cn
    ElementData { outer: 2, valences: &[-1] },         // 113 Nh
    ElementData { outer: 2, valences: &[-1] },         // 114 Fl
    ElementData { outer: 2, valences: &[-1] },         // 115 Mc
    ElementData { outer: 2, valences: &[-1] },         // 116 Lv
    ElementData { outer: 2, valences: &[-1] },         // 117 Ts
    ElementData { outer: 2, valences: &[-1] },         // 118 Og
];
