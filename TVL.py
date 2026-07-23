"""
TVL.py  —  Topological Vortex Logic
=====================================
Classification of stable winding states on the three-torus
T3 = R^3 / L Z^3, with a selected Z3 grading.

This module is deliberately split into two parts:

  MATHEMATICAL CORE  (class TVL, dataclass TVLState)
      Everything here follows from an explicitly specified discrete model: the
      integer winding lattice Z^3, the stipulated quadratic split cost |w|^2,
      and a chosen diagonal homomorphism Z^3 -> Z3.  The stable finite set,
      shell structure, projected weight orbits, and A2/B3 root systems then
      follow exactly.  A TVLState carries only mathematical fields.

  PHYSICAL ADAPTER  (class TVLInterpretation, dataclass PhysicalReading)
      An OPTIONAL, CONJECTURAL reading of the mathematics in Standard-Model
      language: a baryon number and a particle-sector label. These are an
      analogy MATCHED to the mathematics, NOT derived from it. The construction
      contains no spinors, chirality, gauge dynamics, Yukawa couplings, or mass
      spectrum, so nothing in the adapter follows as a theorem. It is provided
      for interpretation only and is kept in a separate object on purpose.

Two distinct order-three structures appear and must not be confused:
    * the CHARGE CLASS   q3 = tr(w) mod 3   -- a linear functional
      H^1(T3,Z) -> Z3 (the chosen "diagonal" functional); see charge_class().
    * the COORDINATE C3 action   g:(w1,w2,w3) -> (w2,w3,w1)   -- a permutation
      of the winding directions, used only for the module invariants; see
      z3_module_invariants().
The target Z3 and the acting cyclic group C3 are distinct structures, even though both have order three.

USAGE
    from TVL import TVL, TVLInterpretation
    s = TVL.classify((1, 0, 0))       # mathematical classification
    print(s)                          #   -> math-only line
    print(s.weight_class)             #   -> 'weight orbit {mu_i}'
    r = TVLInterpretation.read((1,0,0))   # optional conjectural reading
    print(r.sector, r.baryon_number)  #   -> 'quark' 1/3

    python TVL.py                     # run self-test + A2/B3 verification
    python TVL.py 1,0,0 1,1,0         # mathematical classifications only
    python TVL.py --interpret 1,0,0   # add the conjectural reading explicitly
    python TVL.py --all               # all 26 mathematical states
    python TVL.py --all --interpret   # include conjectural readings
    python TVL.py --map               # mathematical map for |w|^2 <= 6

RESULTS exercised by the self-test:
    26 stable states = 3^3 - 1;  shells 6 face / 12 edge / 8 corner;
    charge classes 8 (q3=0) / 9 (q3=1) / 9 (q3=2);  A2 and B3 root systems;
    coordinate-C3 module invariants distinguishing the three shells.

Author: Vladimer Merebashvili
Date:   April 2026 (refactored 2026)
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import permutations, product
from typing import Dict, List, Optional, Tuple

# -- Types ----------------------------------------------------------------------
WindingVector = Tuple[int, int, int]


def _validate(w) -> WindingVector:
    """Return w as a 3-tuple of ints, or raise ValueError.

    Guards against silently mis-classifying inputs of the wrong shape.
    """
    if not isinstance(w, (tuple, list)) or len(w) != 3 or \
            any(type(x) is not int for x in w):
        raise ValueError(
            f"winding vector must be exactly three integers, got {w!r}")
    return (w[0], w[1], w[2])


# ==============================================================================
#  MATHEMATICAL CORE
# ==============================================================================
@dataclass(frozen=True)
class TVLState:
    """Immutable MATHEMATICAL classification of a winding vector w.

    Every field is derived from the explicitly specified mathematical model; there
    is no physical content here. The conjectural Standard-Model reading lives in a
    separate object, TVLInterpretation / PhysicalReading.

    Fields
    ------
    w              : the winding vector (w1, w2, w3) in Z^3
    norm2          : |w|^2 = w1^2 + w2^2 + w3^2
    stable         : True iff |w|^2 <= 3 (Theorem 1)
    shell          : 'vacuum' | 'face' | 'edge' | 'corner' | 'unstable'
    shell_index    : 1 (face) / 2 (edge) / 3 (corner); None otherwise
    q3             : charge class = tr(w) mod 3 in {0,1,2} (chosen functional)
    traceless_norm2: |w_t|^2 for w_t = w - (tr(w)/3)(1,1,1)
    weight_class   : weight-vector content of w_t (NOT an asserted SU(3) irrep)
    """
    w:               WindingVector
    norm2:           int
    stable:          bool
    shell:           str
    shell_index:     Optional[int]
    q3:              int
    traceless_norm2: Fraction
    weight_class:    str

    def __str__(self) -> str:
        if self.shell == 'vacuum':
            return f"w={self.w}  VACUUM"
        if not self.stable:
            return (f"w={str(self.w):>18}  |w|^2={self.norm2}  UNSTABLE "
                    f"(splits favorably; not a stable shell)")
        return (
            f"w={str(self.w):>18}  |w|^2={self.norm2}  "
            f"shell={self.shell:<6}  idx={self.shell_index}  "
            f"q3={self.q3}  |w_t|^2={str(self.traceless_norm2):>4}  "
            f"weight={self.weight_class}"
        )

    def to_dict(self) -> dict:
        """Return the mathematical fields as a plain dictionary."""
        return {
            "w":               self.w,
            "norm2":           self.norm2,
            "stable":          self.stable,
            "shell":           self.shell,
            "shell_index":     self.shell_index,
            "q3":              self.q3,
            "traceless_norm2": str(self.traceless_norm2),
            "weight_class":    self.weight_class,
        }


class TVL:
    """Winding-state classifier with a selected Z3 grading (mathematical core).

    All methods are static; the class is a namespace for the map. Nothing here
    depends on any physical interpretation.

    Layer 1 - STABILITY
        Integer lattice arithmetic on Z^3.
        Theorem 1: w stable  <=>  |w|^2 <= 3.
        Stable direction (proved, no enumeration): for |w|^2 <= 3 each wi is in
        {-1,0,1}, so ui^2 - wi ui >= 0 for every integer ui, hence
        E(u) + E(w-u) - E(w) = 2 sum_i (ui^2 - wi ui) >= 0.
        Unstable direction: |w|^2 >= 4 forces some |wi| >= 2, and the unit split
        s = sign(wi) e_i lowers the cost (w . s = |wi| >= 2 > 1 = |s|^2).

    Layer 2 - CHARGE CLASS
        q3 = tr(w) mod 3, the chosen diagonal functional H^1(T3,Z) -> Z3.
        Distinct from the coordinate C3 permutation action used in Layer 3.

    Layer 3 - WEIGHT CLASSES and ROOT SYSTEMS
        The traceless projection sorts states into weight classes (with q3); it
        does NOT by itself fix an SU(3) irrep. The six traceless edge states
        form the exact A2 root system, and the eighteen face+edge states form
        the type-B3 root system associated with so(7). The three shells are pairwise non-isomorphic
        as rational permutation modules under the coordinate C3 action g.
    """

    # -- Layer 1: stability ----------------------------------------------------

    @staticmethod
    def norm2(w: WindingVector) -> int:
        """Return |w|^2 = w1^2 + w2^2 + w3^2."""
        w = _validate(w)
        return w[0] ** 2 + w[1] ** 2 + w[2] ** 2

    @staticmethod
    def is_stable(w: WindingVector) -> bool:
        """Return True iff w is stable against splitting, i.e. |w|^2 <= 3.

        A split w -> u + (w-u) lowers the stipulated quadratic cost E(w)=|w|^2
        exactly when w . u > |u|^2. See the class docstring for the closed-form
        argument in both directions.
        """
        w = _validate(w)
        return TVL.norm2(w) <= 3

    @staticmethod
    def closed_form_stable(w: WindingVector) -> dict:
        """Run the closed-form stability argument on any winding vector.

        Returns a dict with keys: 'stable', 'norm2', 'reason', 'split_s',
        'w_dot_s'. For |w|^2 >= 4 the returned split_s is an explicit favorable
        unit split; for |w|^2 <= 3 it is None (stable).
        """
        w = _validate(w)
        n2 = TVL.norm2(w)
        if n2 == 0:
            return {'stable': True, 'norm2': 0,
                    'reason': 'Vacuum state. No winding.',
                    'split_s': None, 'w_dot_s': None}
        if n2 <= 3:
            return {'stable': True, 'norm2': n2,
                    'reason': (f'|w|^2={n2} <= 3: stable. Each wi in {{-1,0,1}}, '
                               f'so ui^2 - wi ui >= 0 for every integer ui, '
                               f'hence no split lowers the cost.'),
                    'split_s': None, 'w_dot_s': None}
        for i in range(3):
            if abs(w[i]) >= 2:
                sign_i  = 1 if w[i] > 0 else -1
                s       = tuple(sign_i if j == i else 0 for j in range(3))
                w_dot_s = sum(w[j] * s[j] for j in range(3))
                ws      = tuple(w[j] - s[j] for j in range(3))
                ws_n2   = sum(x ** 2 for x in ws)
                return {'stable': False, 'norm2': n2,
                        'reason': (
                            f'|w|^2={n2} >= 4. Component w[{i}]={w[i]} has '
                            f'|w[{i}]|>=2. Choose s=sign(w[{i}]) e_{i}={s}: '
                            f'w.s={w_dot_s} >= 2 > 1 = |s|^2, so the split is '
                            f'favorable. Residual {ws}, |residual|^2={ws_n2}; '
                            f'cost {n2} -> {1 + ws_n2}.'),
                        'split_s': s, 'w_dot_s': w_dot_s}
        return {'stable': False, 'norm2': n2,
                'reason': 'Unstable (fallback).', 'split_s': None, 'w_dot_s': None}

    @staticmethod
    def shell(w: WindingVector) -> str:
        """Return 'vacuum' / 'face' / 'edge' / 'corner' / 'unstable'."""
        w = _validate(w)
        return {0: 'vacuum', 1: 'face', 2: 'edge', 3: 'corner'}.get(
            TVL.norm2(w), 'unstable')

    @staticmethod
    def shell_index(w: WindingVector) -> Optional[int]:
        """Return the shell index 1 (face) / 2 (edge) / 3 (corner).

        None for vacuum or unstable states. These are the three stable
        norm-shells; stability (|w|^2 <= 3) admits no fourth. (Reading the shell
        index as a Standard-Model generation is a conjectural step handled by
        TVLInterpretation, not a fact of this classifier.)
        """
        w = _validate(w)
        return {1: 1, 2: 2, 3: 3}.get(TVL.norm2(w), None)

    # -- Layer 2: charge class -------------------------------------------------

    @staticmethod
    def charge_class(w: WindingVector) -> int:
        """Return the charge class q3 = (w1 + w2 + w3) mod 3 in {0,1,2}.

        This is a homomorphism H^1(T3,Z) -> Z3 -- specifically the chosen
        "diagonal" functional tr(w) mod 3. A different nonzero functional (a
        unimodular change of basis / shear of T3) would give the SAME global
        8/9/9 split of the 26 nonzero mod-3 vectors but RELABEL individual
        states, so the per-state charge is basis-dependent; only the 8/9/9
        distribution is invariant.

        NOTE: this charge functional targets Z3; it is distinct from the
        coordinate cyclic group C3 generated by g in z3_module_invariants().
        """
        w = _validate(w)
        return (w[0] + w[1] + w[2]) % 3

    # -- Layer 3: traceless projection, weight classes -------------------------

    @staticmethod
    def traceless_projection_norm2(w: WindingVector) -> Fraction:
        """Return |w_t|^2 where w_t = w - (tr(w)/3)(1,1,1).

        Together with q3 this sorts states into weight classes; it does NOT by
        itself fix an SU(3) irrep (many irreps share a projected radius):
            |w_t|^2 = 0    -> weight 0
            |w_t|^2 = 2/3  -> the two conjugate fundamental-weight orbits
            |w_t|^2 = 2    -> A2 root
            |w_t|^2 = 8/3  -> orbit {+-2mu}: 3 weights of the 6 and 3 of 6bar,
                              NOT a single irrep
        """
        w = _validate(w)
        tr  = sum(w)
        w_t = tuple(Fraction(x) - Fraction(tr, 3) for x in w)
        return sum(x ** 2 for x in w_t)

    @staticmethod
    def weight_class(w: WindingVector) -> str:
        """Return the weight class of w -- the weight-vector content of its
        traceless projection, NOT an asserted SU(3) irrep.

        The projected face/edge states realise the fundamental-weight orbits
        {+-mu_i}; q3 distinguishes the two conjugate orbits.
        The |w_t|^2 = 8/3 shell projects to {+-2mu_i}, which is not the
        weight system of a single irreducible representation. Only the A2 root set is an exact, irrep-free
        geometric structure; attaching an irrep NAME is a conjectural step.
        """
        w = _validate(w)
        n2t = TVL.traceless_projection_norm2(w)
        q3  = TVL.charge_class(w)
        if n2t == 0:
            return 'zero weight'
        if n2t == Fraction(2, 3):
            return 'weight orbit {mu_i}' if q3 == 1 else 'weight orbit {-mu_i}'
        if n2t == 2:
            return 'A2 root'
        if n2t == Fraction(8, 3):
            return 'weight orbit {+-2mu_i} (not one irrep)'
        return f'(|w_t|^2={n2t})'

    # -- Full classification ---------------------------------------------------

    @staticmethod
    def classify(w: WindingVector) -> TVLState:
        """Return the mathematical TVLState of a winding vector w.

        Raises ValueError unless w is exactly three integers. The result holds
        mathematical fields only; for the conjectural physical reading use
        TVLInterpretation.read(w).
        """
        w  = _validate(w)
        n2 = TVL.norm2(w)
        n2t = TVL.traceless_projection_norm2(w)

        if n2 == 0:
            return TVLState(w=w, norm2=0, stable=True, shell='vacuum',
                            shell_index=None, q3=0,
                            traceless_norm2=Fraction(0), weight_class='-')
        if n2 >= 4:
            return TVLState(w=w, norm2=n2, stable=False, shell='unstable',
                            shell_index=None, q3=TVL.charge_class(w),
                            traceless_norm2=n2t, weight_class='-')
        return TVLState(
            w=w, norm2=n2, stable=True, shell=TVL.shell(w),
            shell_index=TVL.shell_index(w), q3=TVL.charge_class(w),
            traceless_norm2=n2t, weight_class=TVL.weight_class(w))

    # -- Enumeration helpers ---------------------------------------------------

    @staticmethod
    def all_stable() -> Dict[WindingVector, TVLState]:
        """Return all 26 stable winding states, keyed by winding vector."""
        return {w: TVL.classify(w)
                for w in product(range(-1, 2), repeat=3)
                if 1 <= TVL.norm2(w) <= 3}

    @staticmethod
    def print_all() -> None:
        """Print all 26 mathematical states grouped by shell."""
        states = TVL.all_stable()
        for shell_val, name in [(1, 'FACE'), (2, 'EDGE'), (3, 'CORNER')]:
            print(f'\n-- {name}  |w|^2={shell_val} '
                  + '-' * (44 - len(name)))
            for _, state in sorted(states.items(),
                                   key=lambda x: (x[1].norm2, x[0])):
                if state.norm2 == shell_val:
                    print(state)

    @staticmethod
    def print_full_map() -> None:
        """Print mathematical classifications for non-vacuum |w|^2 <= 6."""
        print(f'{"w":>18}  {"|w|^2":>5}  {"shell":<8}  {"idx":>3}  '
              f'{"q3":>2}  {"weight_class":<42}')
        print('-' * 90)
        for w in sorted(product(range(-2, 3), repeat=3),
                        key=lambda x: (TVL.norm2(x), x)):
            n2 = TVL.norm2(w)
            if n2 == 0 or n2 > 6:
                continue
            state = TVL.classify(w)
            idx = str(state.shell_index) if state.shell_index else '-'
            print(f'{str(w):>18}  {n2:>5}  {state.shell:<8}  {idx:>3}  '
                  f'{state.q3:>2}  {state.weight_class:<42}')

    # -- Module-theoretic invariants (coordinate C3 action) --------------------

    @staticmethod
    def z3_module_invariants() -> dict:
        """Module invariants under the COORDINATE C3 action
        g:(w1,w2,w3) -> (w2,w3,w1) that distinguish the three shells.

        This is a permutation representation on the formal complex vector space
        C[shell] with one basis vector per state. Over a splitting field (C, or
        Q(omega) with omega a primitive cube root of unity) it decomposes as
        rho0^{m0} (+) rho1^{m1} (+) rho2^{m2}, where rho_k is the 1-dim
        character g |-> omega^k. Over Q it is instead trivial^{m0} (+) W^{m1},
        where W is the 2-dim irreducible rational C3-module (the omega and
        omega^2 characters are conjugate and do not split over Q, so m1 = m2).
        Over F3 it fails to be semisimple and F3 has no primitive cube root, so
        neither decomposition holds.

        Over Q, rho0 is the dimension of the invariant subspace and equals
        the number of C3-orbits. The rational character trace is the fixed-basis-point
        count. These are the finer invariants used by the paper. In characteristic 3,
        use invariant-subspace dimension rather than semisimple multiplicity language.

        Returns a dict with 'face', 'edge', 'corner' keys, each a dict with
        'count', 'fixed_pts', 'invariant_dim', 'rho0' (compatibility alias over
        splitting fields), 'rho_nontrivial' (= m1 = m2), and 'uniform'.
        Non-isomorphism: fixed_pts 0 vs 2 separates corner from face and edge;
        invariant-subspace dimension 2 vs 4 separates face from edge.
        """
        g = lambda w: (w[1], w[2], w[0])

        def invariants(family: List[WindingVector]) -> dict:
            n   = len(family)
            fp  = sum(1 for w in family if g(w) == w)
            fp2 = sum(1 for w in family if g(g(w)) == w)
            r0  = (n + fp + fp2) // 3          # trivial-subrep multiplicity
            r_nt = (n - r0) // 2               # multiplicity of each of rho1,rho2
            return {'count': n, 'fixed_pts': fp, 'invariant_dim': r0,
                    'rho0': r0,  # compatibility name over splitting fields
                    'rho_nontrivial': r_nt,
                    'uniform': (fp == 0 and fp2 == 0 and n % 3 == 0)}

        shells = {name: [w for w in product([-1, 0, 1], repeat=3)
                         if sum(x ** 2 for x in w) == d]
                  for name, d in [('face', 1), ('edge', 2), ('corner', 3)]}
        return {name: invariants(fam) for name, fam in shells.items()}

    @staticmethod
    def print_module_invariants() -> None:
        """Print the coordinate-C3 module-invariant table for the three shells."""
        inv = TVL.z3_module_invariants()
        print('\nCoordinate-C3 module invariants (Paper A, Theorem 5)')
        print('-' * 60)
        print(f'  {"shell":<8} {"N":>4} {"fixed":>6} {"inv.dim":>7} '
              f'{"m1=m2":>6} {"uniform":>8}')
        print('  ' + '-' * 56)
        for name in ['face', 'edge', 'corner']:
            d = inv[name]
            print(f'  {name:<8} {d["count"]:>4} {d["fixed_pts"]:>6} '
                  f'{d["invariant_dim"]:>7} {d["rho_nontrivial"]:>6} '
                  f'{"yes" if d["uniform"] else "no":>8}')
        print()
        print('  face vs corner: fixed pts 0 != 2  -> not isomorphic')
        print('  edge vs corner: fixed pts 0 != 2  -> not isomorphic')
        print('  face vs edge:   inv. dim. 2 != 4  -> not isomorphic')

    # -- A2 root-system verification -------------------------------------------

    @staticmethod
    def verify_a2_root_system(verbose: bool = True) -> bool:
        """Verify the six traceless edge states form the A2 root system."""
        roots = [w for w in product((-1, 0, 1), repeat=3)
                 if sum(x * x for x in w) == 2 and sum(w) == 0]
        root_set = set(roots)

        def dot(a, b):
            return sum(a[i] * b[i] for i in range(3))

        def reflect(alpha, beta):
            coefficient = 2 * dot(alpha, beta) // dot(alpha, alpha)
            return tuple(beta[i] - coefficient * alpha[i] for i in range(3))

        negation = all(tuple(-x for x in root) in root_set for root in roots)

        def parallel(a, b):
            return all(a[i] * b[j] == a[j] * b[i]
                       for i in range(3) for j in range(3))

        multiples = all(not parallel(a, b) or a == b
                        or a == tuple(-x for x in b)
                        for a in roots for b in roots)
        integrality = all((2 * dot(a, b)) % dot(a, a) == 0
                          for a in roots for b in roots)
        closure = all(reflect(a, b) in root_set for a in roots for b in roots)
        simple = [(1, -1, 0), (0, 1, -1)]
        cartan = [[2 * dot(a, b) // dot(a, a) for b in simple]
                  for a in simple]
        cartan_ok = cartan == [[2, -1], [-1, 2]]
        ok = (len(roots) == 6 and negation and multiples and integrality
              and closure and cartan_ok)
        if verbose:
            print(f'  Six roots:                              {"pass" if len(roots) == 6 else "FAIL"}')
            print(f'  Negation and reducedness:              {"pass" if negation and multiples else "FAIL"}')
            print(f'  Crystallographic integrality:          {"pass" if integrality else "FAIL"}')
            print(f'  Reflection closure:                    {"pass" if closure else "FAIL"}')
            print(f'  Cartan matrix [[2,-1],[-1,2]]:         {"pass" if cartan_ok else "FAIL"}')
            print(f'  A2: {"ALL CHECKS PASS" if ok else "CHECK FAILURE"}')
        return ok

    # -- B3 root-system verification -------------------------------------------

    @staticmethod
    def verify_b3_root_system(verbose: bool = True) -> bool:
        """Verify that the 18 face-and-edge states form the standard type-B3 root
        system, and that the simple-root reflections generate the order-48
        hyperoctahedral group of signed permutations.

        Checks the four root-system axioms, the B3 Cartan matrix of the simple
        roots (1,-1,0),(0,1,-1),(0,0,1), and -- by explicit generation and
        counting -- that the group of signed 3-coordinate permutations has order
        48 and preserves the root set.
        """
        face = [w for w in product([-1, 0, 1], repeat=3)
                if sum(x ** 2 for x in w) == 1]
        edge = [w for w in product([-1, 0, 1], repeat=3)
                if sum(x ** 2 for x in w) == 2]
        R, R_set = face + edge, set(face + edge)
        ok = True

        def p(msg: str) -> None:
            if verbose:
                print(msg)

        neg = all(tuple(-x for x in w) in R_set for w in R)
        p(f'  Axiom 1 (closure under negation):  {"pass" if neg else "FAIL"}')
        ok = ok and neg

        def parallel(a, b):
            return all(a[i] * b[j] == a[j] * b[i]
                       for i in range(3) for j in range(3))

        mult = all(not parallel(a, b) or a == b
                   or a == tuple(-x for x in b)
                   for a in R for b in R)
        p(f'  Axiom 2 (no non-+-1 multiples):    {"pass" if mult else "FAIL"}')
        ok = ok and mult

        def dot(a, b): return sum(a[i] * b[i] for i in range(3))
        cartan = all((2 * dot(a, b)) % dot(a, a) == 0
                     for a in R for b in R if a != b)
        p(f'  Axiom 3 (integer Cartan entries, {len(R)*(len(R)-1)} pairs):  '
          f'{"pass" if cartan else "FAIL"}')
        ok = ok and cartan

        def reflect(a, b):
            c = 2 * dot(a, b) // dot(a, a)
            return tuple(b[i] - c * a[i] for i in range(3))
        refl = all(reflect(a, b) in R_set for a in R for b in R if a != b)
        p(f'  Axiom 4 (reflection closure):      {"pass" if refl else "FAIL"}')
        ok = ok and refl

        simple   = [(1, -1, 0), (0, 1, -1), (0, 0, 1)]
        A        = [[2 * dot(a, b) // dot(a, a) for b in simple] for a in simple]
        expected = [[2, -1, 0], [-1, 2, -1], [0, -2, 2]]
        cart_ok  = (A == expected)
        p(f'  Cartan matrix of simple roots:     {"pass" if cart_ok else "FAIL"}')
        if verbose and not cart_ok:
            print(f'    computed = {A}\n    expected = {expected}')
        ok = ok and cart_ok

        # Weyl group: generate the matrix group from the three simple-root
        # reflections, then compare it with all signed coordinate permutations.
        def reflection_matrix(alpha):
            columns = []
            a2 = dot(alpha, alpha)
            for j in range(3):
                basis = tuple(1 if i == j else 0 for i in range(3))
                coefficient = 2 * dot(alpha, basis) // a2
                columns.append(tuple(basis[i] - coefficient * alpha[i]
                                     for i in range(3)))
            return tuple(tuple(columns[j][i] for j in range(3))
                         for i in range(3))

        def multiply(A, B):
            return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(3))
                               for j in range(3)) for i in range(3))

        def apply(M, v):
            return tuple(sum(M[i][k] * v[k] for k in range(3)) for i in range(3))

        identity = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        generators = [reflection_matrix(alpha) for alpha in simple]
        generated = {identity}
        frontier = [identity]
        while frontier:
            current = frontier.pop()
            for generator in generators:
                candidate = multiply(generator, current)
                if candidate not in generated:
                    generated.add(candidate)
                    frontier.append(candidate)

        signed_permutations = set()
        for perm in permutations(range(3)):
            for signs in product((1, -1), repeat=3):
                signed_permutations.add(
                    tuple(tuple(signs[i] if perm[i] == j else 0
                                for j in range(3)) for i in range(3)))

        weyl_ok = (len(generated) == 48
                   and generated == signed_permutations
                   and all(apply(M, root) in R_set
                           for M in generated for root in R))
        p(f'  Simple reflections generate order {len(generated)} = 48 and preserve R:  '
          f'{"pass" if weyl_ok else "FAIL"}')
        ok = ok and weyl_ok

        p(f'  Type B3 (associated with so(7)): {"ALL CHECKS PASS" if ok else "CHECK FAILURE"}')
        return ok


# ==============================================================================
#  PHYSICAL ADAPTER  (optional, conjectural -- NOT derived from the mathematics)
# ==============================================================================
@dataclass(frozen=True)
class PhysicalReading:
    """A conjectural Standard-Model reading of a mathematical TVLState.

    These labels are an ANALOGY matched to the mathematics, not consequences of
    it. `baryon_number` rests on one imported QCD identification (Z3 <-> SU(3)
    colour centre); `sector` is an interpretive particle label. Neither is
    derived from the specified lattice model.
    """
    w:             WindingVector
    baryon_number: Fraction
    sector:        str
    note:          str = ("Conjectural reading, matched to the mathematics, "
                          "not derived from it.")


class TVLInterpretation:
    """Optional adapter supplying the conjectural physical reading.

    Kept separate from TVL on purpose: the mathematical classifier stands on its
    own, and none of the labels below follow from it as theorems.
    """

    @staticmethod
    def _require_stable_nonvacuum(w: WindingVector) -> WindingVector:
        """Validate w and require a stable, non-vacuum mathematical state."""
        w = _validate(w)
        if not (1 <= TVL.norm2(w) <= 3):
            raise ValueError(
                'a physical reading is defined only for stable non-vacuum states')
        return w

    @staticmethod
    def baryon_number(w: WindingVector) -> Fraction:
        """Conjectural baryon number: B = +1/3 (q3=1), -1/3 (q3=2), 0 (q3=0).

        Rests on an imported identification of the selected Z3 grading with the SU(3)
        colour centre (q3 with colour triality). The charge class is mathematical;
        choosing the representatives 0,+-1/3 is an imported QCD convention, since
        the grading determines baryon number only modulo integers.
        """
        w = TVLInterpretation._require_stable_nonvacuum(w)
        return {0: Fraction(0), 1: Fraction(1, 3), 2: Fraction(-1, 3)}[
            TVL.charge_class(w)]

    @staticmethod
    def sector(w: WindingVector) -> str:
        """Conjectural particle-sector label matched to the mathematics.

        {+-2mu_i} states -> 'exotic'; q3=0 states ->
        'lepton' (|w_t|^2=0) or 'gluon' (|w_t|^2=2); q3=1 -> 'quark';
        q3=2 -> 'antiquark'. This 6/6/6/6/2 partition is finer than, and must
        not be confused with, the 8/9/9 charge-class counts.
        """
        w = TVLInterpretation._require_stable_nonvacuum(w)
        q3, n2t = TVL.charge_class(w), TVL.traceless_projection_norm2(w)
        if n2t == Fraction(8, 3):
            return 'exotic'
        if q3 == 0:
            return 'lepton' if n2t == 0 else 'gluon'
        return 'quark' if q3 == 1 else 'antiquark'

    @staticmethod
    def print_all() -> None:
        """Print all stable mathematical states with conjectural readings."""
        states = TVL.all_stable()
        for shell_value, name in [(1, 'FACE'), (2, 'EDGE'), (3, 'CORNER')]:
            print(f'\n-- {name}  |w|^2={shell_value} '
                  + '-' * (44 - len(name)))
            for w, state in sorted(states.items(),
                                   key=lambda x: (x[1].norm2, x[0])):
                if state.norm2 == shell_value:
                    reading = TVLInterpretation.read(w)
                    print(f'{state}  [conjectural: {reading.sector}, '
                          f'B={reading.baryon_number}]')

    @staticmethod
    def print_full_map() -> None:
        """Print |w|^2 <= 6, adding readings only for stable states."""
        print(f'{"w":>18}  {"|w|^2":>5}  {"shell":<8}  {"idx":>3}  '
              f'{"q3":>2}  {"weight_class":<42}  [conjectural]')
        print('-' * 108)
        for w in sorted(product(range(-2, 3), repeat=3),
                        key=lambda x: (TVL.norm2(x), x)):
            n2 = TVL.norm2(w)
            if n2 == 0 or n2 > 6:
                continue
            state = TVL.classify(w)
            idx = str(state.shell_index) if state.shell_index else '-'
            reading = TVLInterpretation.read(w)
            label = (f'{reading.sector}, B={reading.baryon_number}'
                     if reading else '-')
            print(f'{str(w):>18}  {n2:>5}  {state.shell:<8}  {idx:>3}  '
                  f'{state.q3:>2}  {state.weight_class:<42}  [{label}]')

    @staticmethod
    def read(w: WindingVector) -> Optional[PhysicalReading]:
        """Return the conjectural PhysicalReading of w, or None if w is not a
        stable non-vacuum state (only stable states get a physical reading)."""
        w = _validate(w)
        if not (1 <= TVL.norm2(w) <= 3):
            return None
        return PhysicalReading(w=w,
                               baryon_number=TVLInterpretation.baryon_number(w),
                               sector=TVLInterpretation.sector(w))


# ==============================================================================
#  Self-test
# ==============================================================================
def _self_test() -> bool:
    """Exercise the classifier and the conjectural adapter against known values.

    These are check() calls, not Python assertions, and they exercise the
    mathematical facts plus a sample of interpretive lookups; they are not a
    proof (the analytic proofs live in the paper).
    """
    passed = [0]
    failures: List[str] = []

    def check(label: str, got, expected) -> None:
        if got == expected:
            passed[0] += 1
        else:
            failures.append(f'{label}: got {got!r}, expected {expected!r}')
            print(f'  FAIL {label}: got {got!r}, expected {expected!r}')

    # -- mathematical core: counts --------------------------------------------
    m = TVL.all_stable()
    check('26 stable states', len(m), 26)
    check('6 face',   sum(1 for s in m.values() if s.shell == 'face'), 6)
    check('12 edge',  sum(1 for s in m.values() if s.shell == 'edge'), 12)
    check('8 corner', sum(1 for s in m.values() if s.shell == 'corner'), 8)
    check('8 q3=0', sum(1 for s in m.values() if s.q3 == 0), 8)
    check('9 q3=1', sum(1 for s in m.values() if s.q3 == 1), 9)
    check('9 q3=2', sum(1 for s in m.values() if s.q3 == 2), 9)
    check('6 A2-root states',
          sum(1 for s in m.values() if s.weight_class == 'A2 root'), 6)
    check('6 exotic-weight states',
          sum(1 for s in m.values()
              if s.weight_class.startswith('weight orbit {+-2mu_i}')), 6)

    # -- mathematical core: representative states ------------------------------
    s = TVL.classify((1, 0, 0))
    check('(1,0,0) stable', s.stable, True)
    check('(1,0,0) shell', s.shell, 'face')
    check('(1,0,0) shell_index', s.shell_index, 1)
    check('(1,0,0) q3', s.q3, 1)
    check('(1,0,0) weight_class', s.weight_class, 'weight orbit {mu_i}')
    check('(1,0,0) traceless_norm2', s.traceless_norm2, Fraction(2, 3))

    s = TVL.classify((1, 1, 1))
    check('(1,1,1) shell', s.shell, 'corner')
    check('(1,1,1) shell_index', s.shell_index, 3)
    check('(1,1,1) q3', s.q3, 0)
    check('(1,1,1) weight_class', s.weight_class, 'zero weight')

    s = TVL.classify((1, -1, 0))
    check('(1,-1,0) q3', s.q3, 0)
    check('(1,-1,0) weight_class', s.weight_class, 'A2 root')
    check('(1,-1,0) traceless_norm2', s.traceless_norm2, Fraction(2))

    for wv in [(1, 1, -1), (-1, -1, 1), (1, -1, 1), (-1, 1, -1)]:
        s = TVL.classify(wv)
        check(f'{wv} weight_class',
              s.weight_class, 'weight orbit {+-2mu_i} (not one irrep)')

    s = TVL.classify((-1, -1, 0))
    check('(-1,-1,0) weight_class', s.weight_class, 'weight orbit {mu_i}')
    s = TVL.classify((1, 1, 0))
    check('(1,1,0) weight_class', s.weight_class, 'weight orbit {-mu_i}')

    check('(2,0,0) unstable', TVL.classify((2, 0, 0)).stable, False)
    check('(2,0,0) shell', TVL.classify((2, 0, 0)).shell, 'unstable')
    check('(0,0,0) vacuum', TVL.classify((0, 0, 0)).shell, 'vacuum')

    # -- input validation on every public vector helper -----------------------
    helpers = [TVL.norm2, TVL.is_stable, TVL.closed_form_stable, TVL.shell,
               TVL.shell_index, TVL.charge_class,
               TVL.traceless_projection_norm2, TVL.weight_class, TVL.classify]
    for helper in helpers:
        rejected = True
        for bad in [(1, 0, 0, 9), (1, 0), (1.0, 0, 0), 'abc']:
            try:
                helper(bad)
                rejected = False
            except ValueError:
                pass
        check(f'{helper.__name__} rejects malformed vectors', rejected, True)

    # -- closed-form stability -------------------------------------------------
    cf = TVL.closed_form_stable((2, 0, 0))
    check('cf(2,0,0) stable', cf['stable'], False)
    check('cf(2,0,0) split_s', cf['split_s'], (1, 0, 0))
    check('cf(2,0,0) w_dot_s', cf['w_dot_s'], 2)
    check('cf(1,1,1) stable', TVL.closed_form_stable((1, 1, 1))['stable'], True)
    check('cf(-3,0,0) w_dot_s', TVL.closed_form_stable((-3, 0, 0))['w_dot_s'], 3)
    check('all 26 states stable by closed form',
          all(TVL.closed_form_stable(w)['stable'] for w in m), True)
    unstable_cube = [w for w in product(range(-3, 4), repeat=3)
                     if TVL.norm2(w) >= 4]
    check('all unstable cube states return a favorable witness',
          all((not TVL.closed_form_stable(w)['stable']
               and TVL.closed_form_stable(w)['split_s'] is not None)
              for w in unstable_cube), True)

    # -- module invariants -----------------------------------------------------
    inv = TVL.z3_module_invariants()
    check('face fixed_pts', inv['face']['fixed_pts'], 0)
    check('corner fixed_pts', inv['corner']['fixed_pts'], 2)
    check('face invariant_dim', inv['face']['invariant_dim'], 2)
    check('edge invariant_dim', inv['edge']['invariant_dim'], 4)
    check('corner invariant_dim', inv['corner']['invariant_dim'], 4)

    # -- conjectural adapter (kept separate) -----------------------------------
    check('interp(1,0,0) sector', TVLInterpretation.read((1, 0, 0)).sector, 'quark')
    check('interp(1,0,0) B', TVLInterpretation.read((1, 0, 0)).baryon_number,
          Fraction(1, 3))
    check('interp(1,1,1) sector', TVLInterpretation.read((1, 1, 1)).sector, 'lepton')
    check('interp(1,-1,0) sector', TVLInterpretation.read((1, -1, 0)).sector, 'gluon')
    check('interp(1,1,-1) sector', TVLInterpretation.read((1, 1, -1)).sector, 'exotic')
    check('interp(1,1,0) sector', TVLInterpretation.read((1, 1, 0)).sector, 'antiquark')
    check('interp(2,0,0) is None', TVLInterpretation.read((2, 0, 0)), None)
    adapter_rejects = True
    for method in [TVLInterpretation.baryon_number, TVLInterpretation.sector]:
        for invalid in [(0, 0, 0), (2, 0, 0)]:
            try:
                method(invalid)
                adapter_rejects = False
            except ValueError:
                pass
    check('direct adapter helpers reject vacuum/unstable states',
          adapter_rejects, True)

    # -- root systems ----------------------------------------------------------
    check('A2 root system', TVL.verify_a2_root_system(verbose=False), True)
    check('B3 root system', TVL.verify_b3_root_system(verbose=False), True)

    if failures:
        print(f'\n{len(failures)} check(s) FAILED.')
        return False
    print(f'All {passed[0]} checks passed '
          f'(mathematical core + conjectural adapter + input validation).')
    return True


# ==============================================================================
#  Entry point
# ==============================================================================
if __name__ == '__main__':
    import sys

    if len(sys.argv) == 1 or '--test' in sys.argv:
        print('Running TVL self-test...\n')
        ok = _self_test()
        print('\n-- A2 Root System Verification ' + '-' * 30)
        a2ok = TVL.verify_a2_root_system()
        print('\n-- B3 Root System Verification ' + '-' * 30)
        b3ok = TVL.verify_b3_root_system()
        print('\n-- Coordinate-C3 Module Invariants ' + '-' * 26)
        TVL.print_module_invariants()
        sys.exit(0 if (ok and a2ok and b3ok) else 1)

    elif '--all' in sys.argv:
        (TVLInterpretation.print_all() if '--interpret' in sys.argv
         else TVL.print_all())

    elif '--map' in sys.argv:
        (TVLInterpretation.print_full_map() if '--interpret' in sys.argv
         else TVL.print_full_map())

    else:
        include_interpretation = '--interpret' in sys.argv
        for arg in sys.argv[1:]:
            if arg.startswith('--'):
                continue
            try:
                cleaned = arg.replace('(', '').replace(')', '')
                parts   = [p for p in cleaned.replace(' ', ',').split(',')
                           if p.strip()]
                w       = tuple(int(p.strip()) for p in parts)
                s = TVL.classify(w)          # raises ValueError unless 3 ints
                reading = (TVLInterpretation.read(w)
                           if include_interpretation else None)
                interpretation = (
                    f'  [conjectural: {reading.sector}, '
                    f'B={reading.baryon_number}]' if reading else '')
                print(f'{s}{interpretation}')
            except (ValueError, TypeError) as e:
                print(f'Could not classify {arg!r}: {e}')
