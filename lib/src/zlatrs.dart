import 'dart:math';

import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dzasum.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/blas/izamax.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zaxpy.dart';
import 'package:lapack/src/blas/zdotc.dart';
import 'package:lapack/src/blas/zdotu.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/ztrsv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zladiv.dart';

void zlatrs(
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final String NORMIN,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> X_,
  final Box<double> SCALE,
  final Array<double> CNORM_,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final X = X_.having();
  final CNORM = CNORM_.having();

  const ZERO = 0.0, HALF = 0.5, ONE = 1.0, TWO = 2.0;
  bool NOTRAN, NOUNIT, UPPER;
  int I, IMAX, J, JFIRST, JINC, JLAST;
  double BIGNUM, GROW, REC, SMLNUM, TJJ, TMAX, TSCAL, XBND, XJ, XMAX;
  Complex CSUMJ, TJJS = Complex.zero, USCAL;

  double CABS1(Complex ZDUM) => ZDUM.toDouble().abs() + ZDUM.imaginary.abs();
  double CABS2(Complex ZDUM) =>
      (ZDUM.toDouble() / 2.0).abs() + (ZDUM.imaginary / 2.0).abs();

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  NOTRAN = lsame(TRANS, 'N');
  NOUNIT = lsame(DIAG, 'N');

  // Test the input parameters.

  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (!NOTRAN && !lsame(TRANS, 'T') && !lsame(TRANS, 'C')) {
    INFO.value = -2;
  } else if (!NOUNIT && !lsame(DIAG, 'U')) {
    INFO.value = -3;
  } else if (!lsame(NORMIN, 'Y') && !lsame(NORMIN, 'N')) {
    INFO.value = -4;
  } else if (N < 0) {
    INFO.value = -5;
  } else if (LDA < max(1, N)) {
    INFO.value = -7;
  }
  if (INFO.value != 0) {
    xerbla('ZLATRS', -INFO.value);
    return;
  }

  // Quick return if possible

  SCALE.value = ONE;
  if (N == 0) return;

  // Determine machine dependent parameters to control overflow.

  SMLNUM = dlamch('Safe minimum') / dlamch('Precision');
  BIGNUM = ONE / SMLNUM;

  if (lsame(NORMIN, 'N')) {
    // Compute the 1-norm of each column, not including the diagonal.

    if (UPPER) {
      // A is upper triangular.

      for (J = 1; J <= N; J++) {
        // 10
        CNORM[J] = dzasum(J - 1, A(1, J).asArray(), 1);
      } // 10
    } else {
      // A is lower triangular.

      for (J = 1; J <= N - 1; J++) {
        // 20
        CNORM[J] = dzasum(N - J, A(J + 1, J).asArray(), 1);
      } // 20
      CNORM[N] = ZERO;
    }
  }

  // Scale the column norms by TSCAL if the maximum element in CNORM is
  // greater than BIGNUM/2.

  IMAX = idamax(N, CNORM, 1);
  TMAX = CNORM[IMAX];
  if (TMAX <= BIGNUM * HALF) {
    TSCAL = ONE;
  } else {
    // Avoid NaN generation if entries in CNORM exceed the
    // overflow threshold

    if (TMAX <= dlamch('Overflow')) {
      // Case 1: All entries in CNORM are valid floating-point numbers
      TSCAL = HALF / (SMLNUM * TMAX);
      dscal(N, TSCAL, CNORM, 1);
    } else {
      // Case 2: At least one column norm of A cannot be
      // represented as a floating-point number. Find the
      // maximum offdiagonal absolute value
      // max( |Re(A(I,J))|, |Im(A(I,J)| ). If this entry is
      // not +/- Infinity, use this value as TSCAL.
      TMAX = ZERO;
      if (UPPER) {
        // A is upper triangular.

        for (J = 2; J <= N; J++) {
          for (I = 1; I <= J - 1; I++) {
            TMAX = max(
                TMAX, max(A[I][J].toDouble().abs(), A[I][J].imaginary.abs()));
          }
        }
      } else {
        // A is lower triangular.

        for (J = 1; J <= N - 1; J++) {
          for (I = J + 1; I <= N; I++) {
            TMAX = max(
                TMAX, max(A[I][J].toDouble().abs(), A[I][J].imaginary.abs()));
          }
        }
      }

      if (TMAX <= dlamch('Overflow')) {
        TSCAL = ONE / (SMLNUM * TMAX);
        for (J = 1; J <= N; J++) {
          if (CNORM[J] <= dlamch('Overflow')) {
            CNORM[J] *= TSCAL;
          } else {
            // Recompute the 1-norm of each column without
            // introducing Infinity in the summation.
            TSCAL = TWO * TSCAL;
            CNORM[J] = ZERO;
            if (UPPER) {
              for (I = 1; I <= J - 1; I++) {
                CNORM[J] += TSCAL * CABS2(A[I][J]);
              }
            } else {
              for (I = J + 1; I <= N; I++) {
                CNORM[J] += TSCAL * CABS2(A[I][J]);
              }
            }
            TSCAL *= HALF;
          }
        }
      } else {
        // At least one entry of A is not a valid floating-point
        // entry. Rely on TRSV to propagate Inf and NaN.
        ztrsv(UPLO, TRANS, DIAG, N, A, LDA, X, 1);
        return;
      }
    }
  }

  // Compute a bound on the computed solution vector to see if the
  // Level 2 BLAS routine ZTRSV can be used.

  XMAX = ZERO;
  for (J = 1; J <= N; J++) {
    // 30
    XMAX = max(XMAX, CABS2(X[J]));
  } // 30
  XBND = XMAX;

  if (NOTRAN) {
    // Compute the growth in A * x = b.

    if (UPPER) {
      JFIRST = N;
      JLAST = 1;
      JINC = -1;
    } else {
      JFIRST = 1;
      JLAST = N;
      JINC = 1;
    }

    if (TSCAL != ONE) {
      GROW = ZERO;
    } else if (NOUNIT) {
      // A is non-unit triangular.

      // Compute GROW = 1/G(j) and XBND = 1/M(j).
      // Initially, G(0) = max{x(i), i=1,...,n}.
      var isTooSmall = false;
      GROW = HALF / max(XBND, SMLNUM);
      XBND = GROW;
      for (J = JFIRST; JINC < 0 ? J >= JLAST : J <= JLAST; J += JINC) {
        // 40

        // Exit the loop if the growth factor is too small.

        if (GROW <= SMLNUM) {
          isTooSmall = true;
          break;
        }

        TJJS = A[J][J];
        TJJ = CABS1(TJJS);

        if (TJJ >= SMLNUM) {
          // M(j) = G(j-1) / abs(A(j,j))

          XBND = min(XBND, min(ONE, TJJ) * GROW);
        } else {
          // M(j) could overflow, set XBND to 0.

          XBND = ZERO;
        }

        if (TJJ + CNORM[J] >= SMLNUM) {
          // G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )

          GROW *= (TJJ / (TJJ + CNORM[J]));
        } else {
          // G(j) could overflow, set GROW to 0.

          GROW = ZERO;
        }
      } // 40
      if (!isTooSmall) {
        GROW = XBND;
      }
    } else {
      // A is unit triangular.

      // Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.

      GROW = min(ONE, HALF / max(XBND, SMLNUM));
      for (J = JFIRST; JINC < 0 ? J >= JLAST : J <= JLAST; J += JINC) {
        // 50

        // Exit the loop if the growth factor is too small.

        if (GROW <= SMLNUM) break;

        // G(j) = G(j-1)*( 1 + CNORM(j) )

        GROW *= (ONE / (ONE + CNORM[J]));
      } // 50
    }
  } else {
    // Compute the growth in A**T * x = b  or  A**H * x = b.

    if (UPPER) {
      JFIRST = 1;
      JLAST = N;
      JINC = 1;
    } else {
      JFIRST = N;
      JLAST = 1;
      JINC = -1;
    }

    if (TSCAL != ONE) {
      GROW = ZERO;
    } else if (NOUNIT) {
      // A is non-unit triangular.

      // Compute GROW = 1/G(j) and XBND = 1/M(j).
      // Initially, M(0) = max{x(i), i=1,...,n}.
      var isTooSmall = false;
      GROW = HALF / max(XBND, SMLNUM);
      XBND = GROW;
      for (J = JFIRST; JINC < 0 ? J >= JLAST : J <= JLAST; J += JINC) {
        // 70

        // Exit the loop if the growth factor is too small.

        if (GROW <= SMLNUM) {
          isTooSmall = true;
          break;
        }

        // G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )

        XJ = ONE + CNORM[J];
        GROW = min(GROW, XBND / XJ);

        TJJS = A[J][J];
        TJJ = CABS1(TJJS);

        if (TJJ >= SMLNUM) {
          // M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))

          if (XJ > TJJ) XBND *= (TJJ / XJ);
        } else {
          // M(j) could overflow, set XBND to 0.

          XBND = ZERO;
        }
      } // 70
      if (!isTooSmall) {
        GROW = min(GROW, XBND);
      }
    } else {
      // A is unit triangular.

      // Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.

      GROW = min(ONE, HALF / max(XBND, SMLNUM));
      for (J = JFIRST; JINC < 0 ? J >= JLAST : J <= JLAST; J += JINC) {
        // 80

        // Exit the loop if the growth factor is too small.

        if (GROW <= SMLNUM) break;

        // G(j) = ( 1 + CNORM(j) )*G(j-1)

        XJ = ONE + CNORM[J];
        GROW /= XJ;
      } // 80
    }
    //  } // 90
  }

  if ((GROW * TSCAL) > SMLNUM) {
    // Use the Level 2 BLAS solve if the reciprocal of the bound on
    // elements of X is not too small.

    ztrsv(UPLO, TRANS, DIAG, N, A, LDA, X, 1);
  } else {
    // Use a Level 1 BLAS solve, scaling intermediate results.

    if (XMAX > BIGNUM * HALF) {
      // Scale X so that its components are less than or equal to
      // BIGNUM in absolute value.

      SCALE.value = (BIGNUM * HALF) / XMAX;
      zdscal(N, SCALE.value, X, 1);
      XMAX = BIGNUM;
    } else {
      XMAX *= TWO;
    }

    if (NOTRAN) {
      // Solve A * x = b

      for (J = JFIRST; JINC < 0 ? J >= JLAST : J <= JLAST; J += JINC) {
        // 120

        // Compute x(j) = b(j) / A(j,j), scaling x if necessary.
        var scale = true;
        XJ = CABS1(X[J]);
        if (NOUNIT) {
          TJJS = A[J][J] * TSCAL.toComplex();
        } else {
          TJJS = TSCAL.toComplex();
          if (TSCAL == ONE) scale = false;
        }
        if (scale) {
          TJJ = CABS1(TJJS);
          if (TJJ > SMLNUM) {
            // abs(A(j,j)) > SMLNUM:

            if (TJJ < ONE) {
              if (XJ > TJJ * BIGNUM) {
                // Scale x by 1/b(j).

                REC = ONE / XJ;
                zdscal(N, REC, X, 1);
                SCALE.value *= REC;
                XMAX *= REC;
              }
            }
            X[J] = zladiv(X[J], TJJS);
            XJ = CABS1(X[J]);
          } else if (TJJ > ZERO) {
            // 0 < abs(A(j,j)) <= SMLNUM:

            if (XJ > TJJ * BIGNUM) {
              // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
              // to avoid overflow when dividing by A(j,j).

              REC = (TJJ * BIGNUM) / XJ;
              if (CNORM[J] > ONE) {
                // Scale by 1/CNORM(j) to avoid overflow when
                // multiplying x(j) times column j.

                REC /= CNORM[J];
              }
              zdscal(N, REC, X, 1);
              SCALE.value *= REC;
              XMAX *= REC;
            }
            X[J] = zladiv(X[J], TJJS);
            XJ = CABS1(X[J]);
          } else {
            // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
            // scale = 0, and compute a solution to A*x = 0.

            for (I = 1; I <= N; I++) {
              // 100
              X[I] = Complex.zero;
            } // 100
            X[J] = Complex.one;
            XJ = ONE;
            SCALE.value = ZERO;
            XMAX = ZERO;
          }
        } // 110

        // Scale x if necessary to avoid overflow when adding a
        // multiple of column j of A.

        if (XJ > ONE) {
          REC = ONE / XJ;
          if (CNORM[J] > (BIGNUM - XMAX) * REC) {
            // Scale x by 1/(2*abs(x(j))).

            REC *= HALF;
            zdscal(N, REC, X, 1);
            SCALE.value *= REC;
          }
        } else if (XJ * CNORM[J] > (BIGNUM - XMAX)) {
          // Scale x by 1/2.

          zdscal(N, HALF, X, 1);
          SCALE.value *= HALF;
        }

        if (UPPER) {
          if (J > 1) {
            // Compute the update
            //    x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)

            zaxpy(J - 1, -X[J] * TSCAL.toComplex(), A(1, J).asArray(), 1, X, 1);
            I = izamax(J - 1, X, 1);
            XMAX = CABS1(X[I]);
          }
        } else {
          if (J < N) {
            // Compute the update
            //    x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)

            zaxpy(N - J, -X[J] * TSCAL.toComplex(), A(J + 1, J).asArray(), 1,
                X(J + 1), 1);
            I = J + izamax(N - J, X(J + 1), 1);
            XMAX = CABS1(X[I]);
          }
        }
      } // 120
    } else if (lsame(TRANS, 'T')) {
      // Solve A**T * x = b

      for (J = JFIRST; JINC < 0 ? J >= JLAST : J <= JLAST; J += JINC) {
        // 170

        // Compute x(j) = b(j) - sum A(k,j)*x(k).
        //                       k<>j

        XJ = CABS1(X[J]);
        USCAL = TSCAL.toComplex();
        REC = ONE / max(XMAX, ONE);
        if (CNORM[J] > (BIGNUM - XJ) * REC) {
          // If x(j) could overflow, scale x by 1/(2*XMAX).

          REC *= HALF;
          if (NOUNIT) {
            TJJS = A[J][J] * TSCAL.toComplex();
          } else {
            TJJS = TSCAL.toComplex();
          }
          TJJ = CABS1(TJJS);
          if (TJJ > ONE) {
            // Divide by A(j,j) when scaling x if A(j,j) > 1.

            REC = min(ONE, REC * TJJ);
            USCAL = zladiv(USCAL, TJJS);
          }
          if (REC < ONE) {
            zdscal(N, REC, X, 1);
            SCALE.value *= REC;
            XMAX *= REC;
          }
        }

        CSUMJ = Complex.zero;
        if (USCAL == Complex.one) {
          // If the scaling needed for A in the dot product is 1,
          // call zdotu to perform the dot product.

          if (UPPER) {
            CSUMJ = zdotu(J - 1, A(1, J).asArray(), 1, X, 1);
          } else if (J < N) {
            CSUMJ = zdotu(N - J, A(J + 1, J).asArray(), 1, X(J + 1), 1);
          }
        } else {
          // Otherwise, use in-line code for the dot product.

          if (UPPER) {
            for (I = 1; I <= J - 1; I++) {
              // 130
              CSUMJ += (A[I][J] * USCAL) * X[I];
            } // 130
          } else if (J < N) {
            for (I = J + 1; I <= N; I++) {
              // 140
              CSUMJ += (A[I][J] * USCAL) * X[I];
            } // 140
          }
        }

        if (USCAL == TSCAL.toComplex()) {
          // Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
          // was not used to scale the dotproduct.
          var scale = true;
          X[J] -= CSUMJ;
          XJ = CABS1(X[J]);
          if (NOUNIT) {
            TJJS = A[J][J] * TSCAL.toComplex();
          } else {
            TJJS = TSCAL.toComplex();
            if (TSCAL == ONE) scale = false;
          }

          if (scale) {
            // Compute x(j) /= A(j,j), scaling if necessary.

            TJJ = CABS1(TJJS);
            if (TJJ > SMLNUM) {
              // abs(A(j,j)) > SMLNUM:

              if (TJJ < ONE) {
                if (XJ > TJJ * BIGNUM) {
                  // Scale X by 1/abs(x(j)).

                  REC = ONE / XJ;
                  zdscal(N, REC, X, 1);
                  SCALE.value *= REC;
                  XMAX *= REC;
                }
              }
              X[J] = zladiv(X[J], TJJS);
            } else if (TJJ > ZERO) {
              // 0 < abs(A(j,j)) <= SMLNUM:

              if (XJ > TJJ * BIGNUM) {
                // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.

                REC = (TJJ * BIGNUM) / XJ;
                zdscal(N, REC, X, 1);
                SCALE.value *= REC;
                XMAX *= REC;
              }
              X[J] = zladiv(X[J], TJJS);
            } else {
              // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
              // scale = 0 and compute a solution to A**T *x = 0.

              for (I = 1; I <= N; I++) {
                // 150
                X[I] = Complex.zero;
              } // 150
              X[J] = Complex.one;
              SCALE.value = ZERO;
              XMAX = ZERO;
            }
          } // 160
        } else {
          // Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
          // product has already been divided by 1/A(j,j).

          X[J] = zladiv(X[J], TJJS) - CSUMJ;
        }
        XMAX = max(XMAX, CABS1(X[J]));
      } // 170
    } else {
      // Solve A**H * x = b

      for (J = JFIRST; JINC < 0 ? J >= JLAST : J <= JLAST; J += JINC) {
        // 220

        // Compute x(j) = b(j) - sum A(k,j)*x(k).
        //                       k<>j

        XJ = CABS1(X[J]);
        USCAL = TSCAL.toComplex();
        REC = ONE / max(XMAX, ONE);
        if (CNORM[J] > (BIGNUM - XJ) * REC) {
          // If x(j) could overflow, scale x by 1/(2*XMAX).

          REC *= HALF;
          if (NOUNIT) {
            TJJS = A[J][J].conjugate() * TSCAL.toComplex();
          } else {
            TJJS = TSCAL.toComplex();
          }
          TJJ = CABS1(TJJS);
          if (TJJ > ONE) {
            // Divide by A(j,j) when scaling x if A(j,j) > 1.

            REC = min(ONE, REC * TJJ);
            USCAL = zladiv(USCAL, TJJS);
          }
          if (REC < ONE) {
            zdscal(N, REC, X, 1);
            SCALE.value *= REC;
            XMAX *= REC;
          }
        }

        CSUMJ = Complex.zero;
        if (USCAL == Complex.one) {
          // If the scaling needed for A in the dot product is 1,
          // call zdotc to perform the dot product.

          if (UPPER) {
            CSUMJ = zdotc(J - 1, A(1, J).asArray(), 1, X, 1);
          } else if (J < N) {
            CSUMJ = zdotc(N - J, A(J + 1, J).asArray(), 1, X(J + 1), 1);
          }
        } else {
          // Otherwise, use in-line code for the dot product.

          if (UPPER) {
            for (I = 1; I <= J - 1; I++) {
              // 180
              CSUMJ += (A[I][J].conjugate() * USCAL) * X[I];
            } // 180
          } else if (J < N) {
            for (I = J + 1; I <= N; I++) {
              // 190
              CSUMJ += (A[I][J].conjugate() * USCAL) * X[I];
            } // 190
          }
        }

        if (USCAL == TSCAL.toComplex()) {
          // Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
          // was not used to scale the dotproduct.
          var scale = true;
          X[J] -= CSUMJ;
          XJ = CABS1(X[J]);
          if (NOUNIT) {
            TJJS = A[J][J].conjugate() * TSCAL.toComplex();
          } else {
            TJJS = TSCAL.toComplex();
            if (TSCAL == ONE) scale = false;
          }

          if (scale) {
            // Compute x(j) /= A(j,j), scaling if necessary.

            TJJ = CABS1(TJJS);
            if (TJJ > SMLNUM) {
              // abs(A(j,j)) > SMLNUM:

              if (TJJ < ONE) {
                if (XJ > TJJ * BIGNUM) {
                  // Scale X by 1/abs(x(j)).

                  REC = ONE / XJ;
                  zdscal(N, REC, X, 1);
                  SCALE.value *= REC;
                  XMAX *= REC;
                }
              }
              X[J] = zladiv(X[J], TJJS);
            } else if (TJJ > ZERO) {
              // 0 < abs(A(j,j)) <= SMLNUM:

              if (XJ > TJJ * BIGNUM) {
                // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.

                REC = (TJJ * BIGNUM) / XJ;
                zdscal(N, REC, X, 1);
                SCALE.value *= REC;
                XMAX *= REC;
              }
              X[J] = zladiv(X[J], TJJS);
            } else {
              // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
              // scale = 0 and compute a solution to A**H *x = 0.

              for (I = 1; I <= N; I++) {
                // 200
                X[I] = Complex.zero;
              } // 200
              X[J] = Complex.one;
              SCALE.value = ZERO;
              XMAX = ZERO;
            }
          } // 210
        } else {
          // Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
          // product has already been divided by 1/A(j,j).

          X[J] = zladiv(X[J], TJJS) - CSUMJ;
        }
        XMAX = max(XMAX, CABS1(X[J]));
      } // 220
    }
    SCALE.value /= TSCAL;
  }

  // Scale the column norms by 1/TSCAL for return.

  if (TSCAL != ONE) {
    dscal(N, ONE / TSCAL, CNORM, 1);
  }
}
