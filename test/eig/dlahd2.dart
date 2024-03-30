import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/lsamen.dart';
import 'package:lapack/src/nio.dart';

void dlahd2(final Nout IOUNIT, final String PATH) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  bool CORZ, SORD;
  String C2;

  // if (IOUNIT <= 0) return;
  SORD = lsame(PATH, 'S') || lsame(PATH, 'D');
  CORZ = lsame(PATH, 'C') || lsame(PATH, 'Z');
  if (!SORD && !CORZ) {
    _print9999(IOUNIT, PATH);
  }
  C2 = PATH.substring(1, 3);

  if (lsamen(2, C2, 'HS')) {
    if (SORD) {
      // Real Non-symmetric Eigenvalue Problem:
      IOUNIT.println('\n ${PATH.a3} -- Real Non-symmetric eigenvalue problem');

      // Matrix types
      _printHSMatrixTypes(IOUNIT, 'pairs ', 'pairs ', 'prs.', 'prs.');

      // Tests performed
      _print9984(IOUNIT, 'orthogonal', 'transpose', '\\');
    } else {
      // Complex Non-symmetric Eigenvalue Problem:
      IOUNIT
          .println('\n ${PATH.a3} -- Complex Non-symmetric eigenvalue problem');

      // Matrix types
      _printHSMatrixTypes(IOUNIT, 'e.vals', 'e.vals', 'e.vs', 'e.vs');

      // Tests performed
      _print9984(IOUNIT, 'unitary', 'conj.transp.', '*');
    }
  } else if (lsamen(2, C2, 'ST')) {
    if (SORD) {
      // Real Symmetric Eigenvalue Problem:

      IOUNIT.println('\n ${PATH.a3} -- Real Symmetric eigenvalue problem');

      // Matrix types
      _printSTMatrixTypes(IOUNIT, 'Symmetric');

      // Tests performed
      IOUNIT.println('\n Tests performed:  See sdrvst.f');
    } else {
      // Complex Hermitian Eigenvalue Problem:

      IOUNIT.println('\n ${PATH.a3} -- Complex Hermitian eigenvalue problem');

      // Matrix types
      _printSTMatrixTypes(IOUNIT, 'Hermitian');

      // Tests performed
      IOUNIT.println('\n Tests performed:  See cdrvst.f');
    }
  } else if (lsamen(2, C2, 'SG')) {
    if (SORD) {
      // Real Symmetric Generalized Eigenvalue Problem:

      IOUNIT.println(
          '\n ${PATH.a3} -- Real Symmetric Generalized eigenvalue problem');

      // Matrix types
      _printSGMatrixTypes(IOUNIT, 'Symmetric');

      // Tests performed

      IOUNIT.println(
          '\n Tests performed:   \n( For each pair (A,B), where A is of the given type \n and B is a random well-conditioned matrix. D is \n diagonal, and Z is orthogonal. )\n 1 = DSYGV, with ITYPE=1 and UPLO=\'U\':  | A Z - B Z D | / ( |A| |Z| n ulp )     \n 2 = DSPGV, with ITYPE=1 and UPLO=\'U\':  | A Z - B Z D | / ( |A| |Z| n ulp )     \n 3 = DSBGV, with ITYPE=1 and UPLO=\'U\':  | A Z - B Z D | / ( |A| |Z| n ulp )     \n 4 = DSYGV, with ITYPE=1 and UPLO=\'L\':  | A Z - B Z D | / ( |A| |Z| n ulp )     \n 5 = DSPGV, with ITYPE=1 and UPLO=\'L\':  | A Z - B Z D | / ( |A| |Z| n ulp )     \n 6 = DSBGV, with ITYPE=1 and UPLO=\'L\':  | A Z - B Z D | / ( |A| |Z| n ulp )     ');
      IOUNIT.println(
          ' 7 = DSYGV, with ITYPE=2 and UPLO=\'U\':  | A B Z - Z D | / ( |A| |Z| n ulp )     \n 8 = DSPGV, with ITYPE=2 and UPLO=\'U\':  | A B Z - Z D | / ( |A| |Z| n ulp )     \n 9 = DSPGV, with ITYPE=2 and UPLO=\'L\':  | A B Z - Z D | / ( |A| |Z| n ulp )     \n10 = DSPGV, with ITYPE=2 and UPLO=\'L\':  | A B Z - Z D | / ( |A| |Z| n ulp )     \n11 = DSYGV, with ITYPE=3 and UPLO=\'U\':  | B A Z - Z D | / ( |A| |Z| n ulp )     \n12 = DSPGV, with ITYPE=3 and UPLO=\'U\':  | B A Z - Z D | / ( |A| |Z| n ulp )     \n13 = DSYGV, with ITYPE=3 and UPLO=\'L\':  | B A Z - Z D | / ( |A| |Z| n ulp )     \n14 = DSPGV, with ITYPE=3 and UPLO=\'L\':  | B A Z - Z D | / ( |A| |Z| n ulp )     ');
    } else {
      // Complex Hermitian Generalized Eigenvalue Problem:

      IOUNIT.println(
          '\n ${PATH.a3} -- Complex Hermitian Generalized eigenvalue problem');

      // Matrix types
      _printSGMatrixTypes(IOUNIT, 'Hermitian');

      // Tests performed
      IOUNIT.println(
          '\n Tests performed:   \n( For each pair (A,B), where A is of the given type \n and B is a random well-conditioned matrix. D is \n diagonal, and Z is unitary. )\n 1 = ZHEGV, with ITYPE=1 and UPLO=\'U\':  | A Z - B Z D | / ( |A| |Z| n ulp )     \n 2 = ZHPGV, with ITYPE=1 and UPLO=\'U\':  | A Z - B Z D | / ( |A| |Z| n ulp )     \n 3 = ZHBGV, with ITYPE=1 and UPLO=\'U\':  | A Z - B Z D | / ( |A| |Z| n ulp )     \n 4 = ZHEGV, with ITYPE=1 and UPLO=\'L\':  | A Z - B Z D | / ( |A| |Z| n ulp )     \n 5 = ZHPGV, with ITYPE=1 and UPLO=\'L\':  | A Z - B Z D | / ( |A| |Z| n ulp )     \n 6 = ZHBGV, with ITYPE=1 and UPLO=\'L\':  | A Z - B Z D | / ( |A| |Z| n ulp )     ');
      IOUNIT.println(
          ' 7 = ZHEGV, with ITYPE=2 and UPLO=\'U\':  | A B Z - Z D | / ( |A| |Z| n ulp )     \n 8 = ZHPGV, with ITYPE=2 and UPLO=\'U\':  | A B Z - Z D | / ( |A| |Z| n ulp )     \n 9 = ZHPGV, with ITYPE=2 and UPLO=\'L\':  | A B Z - Z D | / ( |A| |Z| n ulp )     \n10 = ZHPGV, with ITYPE=2 and UPLO=\'L\':  | A B Z - Z D | / ( |A| |Z| n ulp )     \n11 = ZHEGV, with ITYPE=3 and UPLO=\'U\':  | B A Z - Z D | / ( |A| |Z| n ulp )     \n12 = ZHPGV, with ITYPE=3 and UPLO=\'U\':  | B A Z - Z D | / ( |A| |Z| n ulp )     \n13 = ZHEGV, with ITYPE=3 and UPLO=\'L\':  | B A Z - Z D | / ( |A| |Z| n ulp )     \n14 = ZHPGV, with ITYPE=3 and UPLO=\'L\':  | B A Z - Z D | / ( |A| |Z| n ulp )     ');
    }
  } else if (lsamen(2, C2, 'BD')) {
    if (SORD) {
      // Real Singular Value Decomposition:

      IOUNIT.println('\n ${PATH.a3} -- Real Singular Value Decomposition');

      // Matrix types
      _printBDMatrixTypes(IOUNIT);

      // Tests performed
      _printBDPerformedTests(IOUNIT, 'orthogonal');
    } else {
      // Complex Singular Value Decomposition:

      IOUNIT.println('\n ${PATH.a3} -- Complex Singular Value Decomposition');

      // Matrix types
      _printBDMatrixTypes(IOUNIT);

      // Tests performed
      _printBDPerformedTests(IOUNIT, 'unitary   ');
    }
  } else if (lsamen(2, C2, 'BB')) {
    if (SORD) {
      // Real General Band reduction to bidiagonal form:

      IOUNIT.println('\n ${PATH.a3} -- Real Band reduc. to bidiagonal form');

      // Matrix types
      _printBBMatrixTypes(IOUNIT);

      // Tests performed
      _printBBPerformedTests(IOUNIT, 'orthogonal');
    } else {
      // Complex Band reduction to bidiagonal form:

      IOUNIT.println('\n ${PATH.a3} -- Complex Band reduc. to bidiagonal form');

      // Matrix types
      _printBBMatrixTypes(IOUNIT);

      // Tests performed
      _printBBPerformedTests(IOUNIT, 'unitary   ');
    }
  } else {
    _print9999(IOUNIT, PATH);
    return;
  }
}

void _print9999(final Nout nout, final String path) {
  nout.println(' ${path.a3}:  no header available');
}

void _printHSMatrixTypes(
  final Nout nout,
  final String t1,
  final String t2,
  final String t3,
  final String t4,
) {
  nout.println(' Matrix types (see xCHKHS for details): ');
  nout.println(
      '\n Special Matrices:\n  1=Zero matrix.                          5=Diagonal: geometr. spaced entries.\n  2=Identity matrix.                      6=Diagonal: clustered entries.\n  3=Transposed Jordan block.              7=Diagonal: large, evenly spaced.\n  4=Diagonal: evenly spaced entries.      8=Diagonal: small, evenly spaced.');
  nout.println(
      ' Dense, Non-Symmetric Matrices:\n  9=Well-cond., evenly spaced eigenvals. 14=Ill-cond., geomet. spaced eigenals.\n 10=Well-cond., geom. spaced eigenvals.  15=Ill-conditioned, clustered e.vals.\n 11=Well-conditioned, clustered e.vals.  16=Ill-cond., random complex ${t1.a6}\n 12=Well-cond., random complex ${t2.a6}    17=Ill-cond., large rand. complx ${t3.a4}\n 13=Ill-conditioned, evenly spaced.      18=Ill-cond., small rand. complx ${t4.a4}');
  nout.println(
      ' 19=Matrix with random O(1) entries.     21=Matrix with small random entries.\n 20=Matrix with large random entries.   ');
}

void _print9984(
  final Nout nout,
  final String a,
  final String b,
  final String sym,
) {
  nout.println(
      '\n Tests performed:   (H is Hessenberg, T is Schur, U and Z are $a,\n${' ' * 20}$sym=$b, W is a diagonal matrix of eigenvalues,\n${' ' * 20}L and R are the left and right eigenvector matrices)\n  1 = | A - U H U$sym | / ( |A| n ulp )           2 = | I - U U$sym | / ( n ulp )\n  3 = | H - Z T Z$sym | / ( |H| n ulp )           4 = | I - Z Z$sym | / ( n ulp )\n  5 = | A - UZ T (UZ)$sym | / ( |A| n ulp )       6 = | I - UZ (UZ)$sym | / ( n ulp )\n  7 = | T(e.vects.) - T(no e.vects.) | / ( |T| ulp )\n  8 = | W(e.vects.) - W(no e.vects.) | / ( |W| ulp )\n  9 = | TR - RW | / ( |T| |R| ulp )      10 = | LT - WL | / ( |T| |L| ulp )\n 11= |HX - XW| / (|H| |X| ulp)  (inv.it) 12= |YH - WY| / (|H| |Y| ulp)  (inv.it)');
}

void _printSTMatrixTypes(final Nout nout, final String n) {
  nout.println(' Matrix types (see xDRVST for details): ');
  nout.println(
      '\n Special Matrices:\n  1=Zero matrix.                          5=Diagonal: clustered entries.\n  2=Identity matrix.                      6=Diagonal: large, evenly spaced.\n  3=Diagonal: evenly spaced entries.      7=Diagonal: small, evenly spaced.\n  4=Diagonal: geometr. spaced entries.');
  nout.println(
      ' Dense $n Matrices:\n  8=Evenly spaced eigenvals.             12=Small, evenly spaced eigenvals.\n  9=Geometrically spaced eigenvals.      13=Matrix with random O(1) entries.\n 10=Clustered eigenvalues.               14=Matrix with large random entries.\n 11=Large, evenly spaced eigenvals.      15=Matrix with small random entries.');
}

void _printSGMatrixTypes(final Nout nout, final String t) {
  nout.println(' Matrix types (see xDRVSG for details): ');
  nout.println(
      '\n Special Matrices:\n  1=Zero matrix.                          5=Diagonal: clustered entries.\n  2=Identity matrix.                      6=Diagonal: large, evenly spaced.\n  3=Diagonal: evenly spaced entries.      7=Diagonal: small, evenly spaced.\n  4=Diagonal: geometr. spaced entries.');
  nout.println(
      ' Dense or Banded $t Matrices: \n  8=Evenly spaced eigenvals.          15=Matrix with small random entries.\n  9=Geometrically spaced eigenvals.   16=Evenly spaced eigenvals, KA=1, KB=1.\n 10=Clustered eigenvalues.            17=Evenly spaced eigenvals, KA=2, KB=1.\n 11=Large, evenly spaced eigenvals.   18=Evenly spaced eigenvals, KA=2, KB=2.\n 12=Small, evenly spaced eigenvals.   19=Evenly spaced eigenvals, KA=3, KB=1.\n 13=Matrix with random O(1) entries.  20=Evenly spaced eigenvals, KA=3, KB=2.\n 14=Matrix with large random entries. 21=Evenly spaced eigenvals, KA=3, KB=3.');
}

void _printBDMatrixTypes(final Nout nout) {
  nout.println(
      ' Matrix types (see xCHKBD for details):\n Diagonal matrices:\n   1: Zero${' ' * 28} 5: Clustered entries\n   2: Identity${' ' * 24} 6: Large, evenly spaced entries\n   3: Evenly spaced entries${' ' * 11} 7: Small, evenly spaced entries\n   4: Geometrically spaced entries\n General matrices:\n   8: Evenly spaced sing. vals.${' ' * 7}12: Small, evenly spaced sing vals\n   9: Geometrically spaced sing vals  13: Random, O(1) entries\n  10: Clustered sing. vals.${' ' * 11}14: Random, scaled near overflow\n  11: Large, evenly spaced sing vals  15: Random, scaled near underflow');
}

void _printBDPerformedTests(final Nout nout, final String t) {
  nout.println(
      '\n Test ratios:  (B: bidiagonal, S: diagonal, Q, P, U, and V: ${t.a10}\n${' ' * 16}X: m x nrhs, Y = Q\' X, and Z = U\' Y)');
  nout.println(
      '   1: norm( A - Q B P\' ) / ( norm(A) max(m,n) ulp )\n   2: norm( I - Q\' Q )   / ( m ulp )\n   3: norm( I - P\' P )   / ( n ulp )\n   4: norm( B - U S V\' ) / ( norm(B) min(m,n) ulp )\n   5: norm( Y - U Z )    / ( norm(Z) max(min(m,n),k) ulp )\n   6: norm( I - U\' U )   / ( min(m,n) ulp )\n   7: norm( I - V\' V )   / ( min(m,n) ulp )\n   8: Test ordering of S  (0 if nondecreasing, 1/ulp  otherwise)\n   9: norm( S - S1 )     / ( norm(S) ulp ), where S1 is computed\n${' ' * 43} without computing U and V\'\n  10: Sturm sequence test (0 if sing. vals of B within THRESH of S)\n  11: norm( A - (QU) S (V\' P\') ) / ( norm(A) max(m,n) ulp )\n  12: norm( X - (QU) Z )         / ( |X| max(M,k) ulp )\n  13: norm( I - (QU)\'(QU) )      / ( M ulp )\n  14: norm( I - (V\' P\') (P V) )  / ( N ulp )\n  15: norm( B - U S V\' ) / ( norm(B) min(m,n) ulp )\n  16: norm( I - U\' U )   / ( min(m,n) ulp )\n  17: norm( I - V\' V )   / ( min(m,n) ulp )\n  18: Test ordering of S  (0 if nondecreasing, 1/ulp  otherwise)\n  19: norm( S - S1 )     / ( norm(S) ulp ), where S1 is computed\n${' ' * 43} without computing U and V\'\n  20: norm( B - U S V\' )  / ( norm(B) min(m,n) ulp )  DBDSVX(V,A)\n  21: norm( I - U\' U )    / ( min(m,n) ulp )\n  22: norm( I - V\' V )    / ( min(m,n) ulp )\n  23: Test ordering of S  (0 if nondecreasing, 1/ulp  otherwise)\n  24: norm( S - S1 )      / ( norm(S) ulp ), where S1 is computed\n${' ' * 44} without computing U and V\'\n  25: norm( S - U\' B V ) / ( norm(B) n ulp )  DBDSVX(V,I)\n  26: norm( I - U\' U )    / ( min(m,n) ulp )\n  27: norm( I - V\' V )    / ( min(m,n) ulp )\n  28: Test ordering of S  (0 if nondecreasing, 1/ulp  otherwise)\n  29: norm( S - S1 )      / ( norm(S) ulp ), where S1 is computed\n${' ' * 44} without computing U and V\'\n  30: norm( S - U\' B V ) / ( norm(B) n ulp )  DBDSVX(V,V)\n  31: norm( I - U\' U )    / ( min(m,n) ulp )\n  32: norm( I - V\' V )    / ( min(m,n) ulp )\n  33: Test ordering of S  (0 if nondecreasing, 1/ulp  otherwise)\n  34: norm( S - S1 )      / ( norm(S) ulp ), where S1 is computed\n${' ' * 44} without computing U and V\'');
}

void _printBBMatrixTypes(final Nout nout) {
  nout.println(
      ' Matrix types (see xCHKBB for details):\n Diagonal matrices:\n   1: Zero${' ' * 28} 5: Clustered entries\n   2: Identity${' ' * 24} 6: Large, evenly spaced entries\n   3: Evenly spaced entries${' ' * 11} 7: Small, evenly spaced entries\n   4: Geometrically spaced entries\n General matrices:\n   8: Evenly spaced sing. vals.${' ' * 7}12: Small, evenly spaced sing vals\n   9: Geometrically spaced sing vals  13: Random, O(1) entries\n  10: Clustered sing. vals.${' ' * 11}14: Random, scaled near overflow\n  11: Large, evenly spaced sing vals  15: Random, scaled near underflow');
}

void _printBBPerformedTests(final Nout nout, final String t) {
  nout.println(
      '\n Test ratios:  (B: upper bidiagonal, Q and P: ${t.a10}\n${' ' * 16}C: m x nrhs, PT = P\', Y = Q\' C)\n 1: norm( A - Q B PT ) / ( norm(A) max(m,n) ulp )\n 2: norm( I - Q\' Q )   / ( m ulp )\n 3: norm( I - PT PT\' )   / ( n ulp )\n 4: norm( Y - Q\' C )   / ( norm(Y) max(m,nrhs) ulp )');
}
