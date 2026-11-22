// harvest.m: compute raw Hecke spectra with Magma, to make the site's data stores larger.
//
// Magma does the computer algebra; _scripts/magma_import.py checks what it wrote and
// merges it into the JSON files the pages read. Nothing here writes JSON.
//
//   store      what it computes                                   feeds
//   planar     T_2 on S_2(N, chi) for every even chi mod N, N prime   limits/planar.json
//   quadratic  T_2 on S_2(N, (./N)), N prime = 5 mod 8                limits/quadratic.json
//   weight     T_l on S_k(SL_2(Z)), l = 2, 3, 5                      limits/weight.json
//   spectra    Brandt matrices B(2), B(3), B(5) of G_p(l)             isogeny/spectra.json
//   diagonal   Brandt matrices B(2), B(l), l largest prime <= (log p)^2   limits/diagonal.json
//
// Run from the repository root, for example
//
//   magma -b store:=planar,quadratic out:=harvest.txt _scripts/magma/harvest.m
//   magma -b store:=planar planarN:=769,1009,1201 out:=planar.txt _scripts/magma/harvest.m
//   python3 _scripts/magma_import.py harvest.txt planar.txt
//
// Command-line variables (all optional; Magma passes them as strings):
//   store      comma-separated stores (default: all five)
//   out        output file, appended to (default harvest.txt)
//   planarN    primes N for planar          quadN     primes N = 5 mod 8 for quadratic
//   weightK    weights k for weight         weightL   primes l for weight (default 2,3,5)
//   spectraP   primes p for spectra         diagonalP primes p for diagonal
//   prec       decimal digits for complex roots (default 40)
// Lists are comma-separated integers. The defaults below continue the current data.
//
// Output: one line per record, fields separated by spaces; lines starting with # are
// comments. Numbers are printed with 12 significant digits.
//   planar N o k d Re(chi(2)) Im(chi(2)) Re(l_1) Im(l_1) ... Re(l_d) Im(l_d)
//       one line per even character chi = (orbit representative of order o)^k-conjugate:
//       the eigenvalues l_i of T_2 on S_2(N, chi), not normalised
//   quadratic N d y_1 ... y_d          y = lambda/(i sqrt 2)
//   weight l k d x_1 ... x_d           x = lambda / l^((k-1)/2)
//   brandt p l n i j,v j,v ...         row i of the Hecke operator T_l on
//       SupersingularModule(p), as Magma returns it (acting on the right): its nonzero
//       entries v = T_l[i, j]; one line per row, n rows
//
// Method notes.
//  * Characters: the even characters of order o (o | (N-1)/2) form one Galois orbit.
//    For each o one representative eps with values in Q(zeta_o) is used; the Hecke
//    polynomial of T_2 on CuspidalSubspace(ModularSymbols(eps, 2, +1)) is computed over
//    Q(zeta_o) and mapped into C by zeta_o -> exp(2 pi i k/o) for each k prime to o, which
//    gives the spectrum for each character in the orbit. Its roots are found numerically.
//  * HeckePolynomial is used rather than the characteristic polynomial of HeckeOperator:
//    for prime n the handbook notes it is often much faster (it uses the Deligne bound).
//  * Repeated eigenvalues: polynomials go through SquarefreeFactorization before roots
//    are taken, so multiplicities are exact.
//  * Brandt matrices come from SupersingularModule(p), which uses the Mestre-Oesterle
//    method of graphs; its basis is the supersingular j-invariants.

if not assigned store then store := "planar,quadratic,weight,spectra,diagonal"; end if;
if not assigned out then out := "harvest.txt"; end if;
if not assigned prec then prec := "40"; end if;

// Plain helpers (no closures, no Join, no globals inside functions: all of that differs
// between Magma versions).
function JoinS(seq, sep)
    if #seq eq 0 then return ""; end if;
    if #seq eq 1 then return seq[1]; end if;
    return seq[1] cat &cat[sep cat seq[i] : i in [2..#seq]];
end function;

function IntList(s)
    return [StringToInteger(t) : t in Split(s, ", ")];
end function;

// Primes near a geometric progression from lo to hi, count of them, each = r mod m.
function GeometricPrimes(lo, hi, count, r, m)
    R := RealField(30);
    out := [];
    for t in [0..count - 1] do
        x := Round(Exp(Log(R!lo) + t * (Log(R!hi) - Log(R!lo)) / (count - 1)));
        q := NextPrime(x - 1);
        while q mod m ne r do q := NextPrime(q); end while;
        if #out eq 0 or q gt out[#out] then Append(~out, q); end if;
    end for;
    return out;
end function;

// The defaults: they continue the current data (see the lists in _scripts/limit_laws.py
// and _scripts/isogeny_graphs.py).
if not assigned planarN then
    planarN := "701,809,907,1009,1103,1201,1301,1409,1601,1801,2003";
end if;
if not assigned quadN then
    quadN := JoinS([IntegerToString(q) : q in GeometricPrimes(4000, 20000, 12, 5, 8)], ",");
end if;
if not assigned weightK then
    weightK := JoinS([IntegerToString(k) : k in [984..2400 by 48]], ",");
end if;
if not assigned weightL then weightL := "2,3,5"; end if;
if not assigned spectraP then
    spectraP := JoinS([IntegerToString(q) : q in GeometricPrimes(50000, 100000, 12, 1, 2)], ",");
end if;
if not assigned diagonalP then
    diagonalP := JoinS([IntegerToString(q) : q in GeometricPrimes(20000, 100000, 16, 1, 2)], ",");
end if;

PREC := StringToInteger(prec);
F := Open(out, "a");

function Num(x)
    return Sprint(RealField(12)!x);
end function;

procedure Emit(F, line)
    fprintf F, "%o\n", line;
    Flush(F);
end procedure;

// An element a of Q or of Q(zeta_o) = CyclotomicField(o) in C, by zeta_o -> z (cyc true;
// Eltseq gives the coordinates in the power basis of zeta_o), or by coercion (cyc false).
function EmbedC(a, C, z, cyc)
    if not cyc then return C ! a; end if;
    c := Eltseq(a);
    return &+[C | (C ! c[i]) * z^(i - 1) : i in [1..#c]];
end function;

// The roots in C (a ComplexField or RealField), with multiplicity, of the image of P under
// EmbedC. The squarefree factorisation is exact, so multiplicities are right.
function RootsEmbedded(P, C, z, cyc)
    rts := [];
    for fm in SquarefreeFactorization(P) do
        f := PolynomialRing(C) ! [EmbedC(a, C, z, cyc) : a in Coefficients(fm[1])];
        for r in Roots(f) do
            rts cat:= [r[1] : j in [1..r[2] * fm[2]]];
        end for;
    end for;
    return rts;
end function;

// ---- planar -----------------------------------------------------------------------

procedure Planar(N, F, PREC)
    t0 := Cputime();
    C := ComplexField(PREC);
    R := RealField(PREC);
    h := (N - 1) div 2;
    total := 0;
    for o in Divisors(h) do
        K := o le 2 select RationalField() else CyclotomicField(o);
        G := DirichletGroup(N, K);
        reps := [e : e in Elements(G) | Order(e) eq o and e(N - 1) eq 1];
        error if #reps eq 0, "no even character of order", o, "mod", N;
        eps := reps[1];
        S := CuspidalSubspace(ModularSymbols(eps, 2, +1));
        d := Dimension(S);
        if d eq 0 then continue; end if;
        P := HeckePolynomial(S, 2);
        for k in [k : k in [1..Max(o, 1)] | Gcd(k, o) eq 1] do
            th := 2 * Pi(R) * k / o;
            z := Cos(th) + C.1 * Sin(th);
            cyc := o gt 2;
            lam := RootsEmbedded(P, C, z, cyc);
            error if #lam ne d, "planar: wrong number of roots at N =", N, "o =", o;
            c := EmbedC(eps(2), C, z, cyc);
            fields := [IntegerToString(N), IntegerToString(o), IntegerToString(k), IntegerToString(d),
                       Num(Real(c)), Num(Imaginary(c))];
            for x in lam do fields cat:= [Num(Real(x)), Num(Imaginary(x))]; end for;
            Emit(F, "planar " cat JoinS(fields, " "));
            total +:= d;
        end for;
    end for;
    error if total ne (N - 5) * (N - 7) div 24, "planar: dimensions do not add up to the genus at N =", N;
    printf "planar N = %o: %o eigenvalues, %o s\n", N, total, Cputime(t0);
end procedure;

// ---- quadratic ----------------------------------------------------------------------

procedure Quadratic(N, F, PREC)
    t0 := Cputime();
    error if N mod 8 ne 5 or not IsPrime(N), "quadratic: N must be a prime = 5 mod 8, not", N;
    eps := [e : e in Elements(DirichletGroup(N)) | Order(e) eq 2][1];
    error if eps(2) ne -1, "quadratic: chi(2) is not -1 at N =", N;
    S := CuspidalSubspace(ModularSymbols(eps, 2, +1));
    d := Dimension(S);
    P := HeckePolynomial(S, 2);               // over Q; its roots are +-i a, a real
    c := Coefficients(P);                     // c[i] is the coefficient of x^(i-1)
    m0 := Min([i - 1 : i in [1..#c] | c[i] ne 0]);
    error if exists{i : i in [1..#c] | IsOdd(i - 1 - m0) and c[i] ne 0},
        "quadratic: the Hecke polynomial is not x^m Q(x^2) at N =", N;
    Rx<u> := PolynomialRing(RationalField());
    Q := Rx ! [c[m0 + 1 + 2 * t] : t in [0..(#c - 1 - m0) div 2]];   // roots u = lambda^2 = -a^2
    R := RealField(PREC);
    us := RootsEmbedded(Q, R, R ! 1, false);
    ys := [R | 0 : j in [1..m0]];
    for x in us do
        error if x gt 10^(-PREC div 2), "quadratic: a positive root of Q at N =", N;
        y := Sqrt(Max(-x, 0) / 2);           // lambda = +-i a, y = a / sqrt 2
        ys cat:= [y, -y];
    end for;
    error if #ys ne d, "quadratic: wrong number of roots at N =", N;
    Emit(F, "quadratic " cat JoinS([IntegerToString(N), IntegerToString(d)] cat [Num(y) : y in ys], " "));
    printf "quadratic N = %o: dim %o, %o s\n", N, d, Cputime(t0);
end procedure;

// ---- weight -------------------------------------------------------------------------

procedure Weight(l, k, F, PREC)
    t0 := Cputime();
    S := CuspidalSubspace(ModularSymbols(1, k, +1));
    d := Dimension(S);
    if d eq 0 then return; end if;
    P := HeckePolynomial(S, l);
    digits := Ceiling((k - 1) / 2 * Log(RealField(20) ! l) / Log(RealField(20) ! 10)) + PREC;
    R := RealField(digits);
    scale := Sqrt(R ! l)^(k - 1);
    xs := [r / scale : r in RootsEmbedded(P, R, R ! 1, false)];
    error if #xs ne d, "weight: wrong number of real roots at l =", l, "k =", k;
    Emit(F, "weight " cat JoinS([IntegerToString(l), IntegerToString(k), IntegerToString(d)] cat [Num(x) : x in xs], " "));
    printf "weight l = %o, k = %o: dim %o, %o s\n", l, k, d, Cputime(t0);
end procedure;

// ---- Brandt matrices ------------------------------------------------------------------

function EllFor(p)
    // the largest prime <= (log p)^2, as in _scripts/limit_laws.py
    n := Floor(Log(RealField(30) ! p)^2);
    return IsPrime(n) select n else PreviousPrime(n);
end function;

procedure Brandt(p, ells, F)
    t0 := Cputime();
    M := SupersingularModule(p);
    n := Dimension(M);
    for l in ells do
        T := HeckeOperator(M, l);
        head := JoinS(["brandt", IntegerToString(p), IntegerToString(l), IntegerToString(n)], " ");
        for i in [1..n] do
            row := T[i];
            try
                js := Sort([j : j in Support(row)]);      // a set or a sequence, by version
            catch err
                js := [j : j in [1..n] | row[j] ne 0];
            end try;
            Emit(F, head cat " " cat IntegerToString(i) cat " "
                    cat JoinS([Sprintf("%o,%o", j, row[j]) : j in js], " "));
        end for;
    end for;
    printf "brandt p = %o, l in %o: dim %o, %o s\n", p, ells, n, Cputime(t0);
end procedure;

// ---- main -------------------------------------------------------------------------------

v1, v2, v3 := GetVersion();
Emit(F, Sprintf("# harvest.m, Magma %o.%o-%o, stores %o", v1, v2, v3, store));
stores := Split(store, ", ");

if "planar" in stores then
    for N in IntList(planarN) do Planar(N, F, PREC); end for;
end if;
if "quadratic" in stores then
    for N in IntList(quadN) do Quadratic(N, F, PREC); end for;
end if;
if "weight" in stores then
    for l in IntList(weightL) do
        for k in IntList(weightK) do Weight(l, k, F, PREC); end for;
    end for;
end if;
// Brandt matrices: one SupersingularModule per prime, with every l that is wanted.
want := AssociativeArray();
if "spectra" in stores then
    for p in IntList(spectraP) do want[p] := {2, 3, 5}; end for;
end if;
if "diagonal" in stores then
    for p in IntList(diagonalP) do
        want[p] := (IsDefined(want, p) select want[p] else {}) join {2, EllFor(p)};
    end for;
end if;
for p in Sort(SetToSequence(Keys(want))) do
    Brandt(p, Sort(SetToSequence(want[p])), F);
end for;

Emit(F, "# done");
quit;
