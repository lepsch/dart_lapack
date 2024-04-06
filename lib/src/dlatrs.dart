import 'package:lapack/lapack.dart';

void dlatrs(
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final String NORMIN,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> X_,
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
  const ZERO = 0.0, HALF = 0.5, ONE = 1.0;
  int I, IMAX, J, JFIRST, JINC, JLAST;
  double GROW, REC, SUMJ, TJJ, TJJS = 0, TMAX, TSCAL, USCAL, XBND, XJ, XMAX;
  final WORK = Array<double>(1);

  INFO.value = 0;
  final UPPER = lsame(UPLO, 'U');
  final NOTRAN = lsame(TRANS, 'N');
  final NOUNIT = lsame(DIAG, 'N');

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
    xerbla('DLATRS', -INFO.value);
    return;
  }

  // Quick return if possible
  SCALE.value = ONE;
  if (N == 0) return;


  // Determine machine dependent parameters to control overflow.
  final SMLNUM = dlamch('Safe minimum') / dlamch('Precision');
  final BIGNUM = ONE / SMLNUM;

  if (lsame(NORMIN, 'N')) {
    // Compute the 1-norm of each column, not including the diagonal.

    if (UPPER) {
      // A is upper triangular.
      for (var J = 1; J <= N; J++) {
        CNORM[J] = dasum(J - 1, A(1, J).asArray(), 1);
      }
    } else {
      // A is lower triangular.
      for (var J = 1; J <= N - 1; J++) {
        CNORM[J] = dasum(N - J, A(J + 1, J).asArray(), 1);
      }
      CNORM[N] = ZERO;
    }
  }

  // Scale the column norms by TSCAL if the maximum element in CNORM is
  // greater than BIGNUM.
  IMAX = idamax(N, CNORM, 1);
  TMAX = CNORM[IMAX];
  if (TMAX <= BIGNUM) {
    TSCAL = ONE;
  } else {
    // Avoid NaN generation if entries in CNORM exceed the
    // overflow threshold

    if (TMAX <= dlamch('Overflow')) {
      // Case 1: All entries in CNORM are valid floating-point numbers
      TSCAL = ONE / (SMLNUM * TMAX);
      dscal(N, TSCAL, CNORM, 1);
    } else {
      // Case 2: At least one column norm of A cannot be represented
      // as floating-point number. Find the offdiagonal entry A[ I][ J ]
      // with the largest absolute value. If this entry is not +/- Infinity,
      // use this value as TSCAL.
      TMAX = ZERO;
      if (UPPER) {
        // A is upper triangular.
        for (var J = 2; J <= N; J++) {
          TMAX = max(dlange('M', J - 1, 1, A(1, J), 1, WORK), TMAX);
        }
      } else {
        // A is lower triangular.
        for (var J = 1; J <= N - 1; J++) {
          TMAX = max(dlange('M', N - J, 1, A(J + 1, J), 1, WORK), TMAX);
        }
      }

      if (TMAX <= dlamch('Overflow')) {
        TSCAL = ONE / (SMLNUM * TMAX);
        for (var J = 1; J <= N; J++) {
          if (CNORM[J] <= dlamch('Overflow')) {
            CNORM[J] *= TSCAL;
          } else {
            // Recompute the 1-norm without introducing Infinity
            // in the summation
            CNORM[J] = ZERO;
            if (UPPER) {
              for (var I = 1; I <= J - 1; I++) {
                CNORM[J] += TSCAL * A[I][J].abs();
              }
            } else {
              for (var I = J + 1; I <= N; I++) {
                CNORM[J] += TSCAL * A[I][J].abs();
              }
            }
          }
        }
      } else {
        // At least one entry of A is not a valid floating-point entry.
        // Rely on TRSV to propagate Inf and NaN.
        dtrsv(UPLO, TRANS, DIAG, N, A, LDA, X, 1);
        return;
      }
    }
  }

  // Compute a bound on the computed solution vector to see if the
  // Level 2 BLAS routine DTRSV can be used.

  J = idamax(N, X, 1);
  XMAX = X[J].abs();
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
      // Initially, G(0) = max{x[i], i=1,...,n}.

      GROW = ONE / max(XBND, SMLNUM);
      XBND = GROW;
      var isTooSmall = false;
      for (var J = JFIRST; JINC < 0 ? J >= JLAST : J <= JLAST; J += JINC) {
        // Exit the loop if the growth factor is too small.

        if (GROW <= SMLNUM) {
          isTooSmall = true;
          break;
        }

        // M(j) = G(j-1) / abs(A[j][j])

        TJJ = A[J][J].abs();
        XBND = min(XBND, min(ONE, TJJ) * GROW);
        if (TJJ + CNORM[J] >= SMLNUM) {
          // G(j) = G(j-1)*( 1 + CNORM[j] / abs(A[j][j]) )

          GROW *= (TJJ / (TJJ + CNORM[J]));
        } else {
          // G(j) could overflow, set GROW to 0.

          GROW = ZERO;
        }
      }
      if (!isTooSmall) {
        GROW = XBND;
      }
    } else {
      // A is unit triangular.

      // Compute GROW = 1/G(j), where G(0) = max{x[i], i=1,...,n}.

      GROW = min(ONE, ONE / max(XBND, SMLNUM));
      for (var J = JFIRST; JINC < 0 ? J >= JLAST : J <= JLAST; J += JINC) {
        // Exit the loop if the growth factor is too small.

        if (GROW <= SMLNUM) break;

        // G(j) = G(j-1)*( 1 + CNORM[j] )

        GROW *= (ONE / (ONE + CNORM[J]));
      }
    }
  } else {
    // Compute the growth in A**T * x = b.

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
      // Initially, M(0) = max{x[i], i=1,...,n}.

      GROW = ONE / max(XBND, SMLNUM);
      XBND = GROW;
      var isTooSmall = false;
      for (var J = JFIRST; JINC < 0 ? J >= JLAST : J <= JLAST; J += JINC) {
        // Exit the loop if the growth factor is too small.

        if (GROW <= SMLNUM) {
          isTooSmall = true;
          break;
        }

        // G(j) = max( G(j-1), M(j-1)*( 1 + CNORM[j] ) )

        XJ = ONE + CNORM[J];
        GROW = min(GROW, XBND / XJ);

        // M(j) = M(j-1)*( 1 + CNORM[j] ) / abs(A[j][j])

        TJJ = A[J][J].abs();
        if (XJ > TJJ) XBND *= (TJJ / XJ);
      }
      if (!isTooSmall) {
        GROW = min(GROW, XBND);
      }
    } else {
      // A is unit triangular.

      // Compute GROW = 1/G(j), where G(0) = max{x[i], i=1,...,n}.

      GROW = min(ONE, ONE / max(XBND, SMLNUM));
      for (var J = JFIRST; JINC < 0 ? J >= JLAST : J <= JLAST; J += JINC) {
        // Exit the loop if the growth factor is too small.

        if (GROW <= SMLNUM) break;

        // G(j) = ( 1 + CNORM[j] )*G(j-1)

        XJ = ONE + CNORM[J];
        GROW /= XJ;
      }
    }
  }

  if ((GROW * TSCAL) > SMLNUM) {
    // Use the Level 2 BLAS solve if the reciprocal of the bound on
    // elements of X is not too small.

    dtrsv(UPLO, TRANS, DIAG, N, A, LDA, X, 1);
  } else {
    // Use a Level 1 BLAS solve, scaling intermediate results.

    if (XMAX > BIGNUM) {
      // Scale X so that its components are less than or equal to
      // BIGNUM in absolute value.

      SCALE.value = BIGNUM / XMAX;
      dscal(N, SCALE.value, X, 1);
      XMAX = BIGNUM;
    }

    if (NOTRAN) {
      // Solve A * x = b

      for (var J = JFIRST; JINC < 0 ? J >= JLAST : J <= JLAST; J += JINC) {
        // Compute x[j] = b(j) / A[j][j], scaling x if necessary.

        XJ = X[J].abs();
        var scale = true;
        if (NOUNIT) {
          TJJS = A[J][J] * TSCAL;
        } else {
          TJJS = TSCAL;
          if (TSCAL == ONE) scale = false;
        }

        if (scale) {
          TJJ = TJJS.abs();
          if (TJJ > SMLNUM) {
            // abs(A[j][j]) > SMLNUM:

            if (TJJ < ONE) {
              if (XJ > TJJ * BIGNUM) {
                // Scale x by 1/b(j).

                REC = ONE / XJ;
                dscal(N, REC, X, 1);
                SCALE.value *= REC;
                XMAX *= REC;
              }
            }
            X[J] /= TJJS;
            XJ = X[J].abs();
          } else if (TJJ > ZERO) {
            // 0 < abs(A[j][j]) <= SMLNUM:

            if (XJ > TJJ * BIGNUM) {
              // Scale x by (1/abs(x[j]))*abs(A[j][j])*BIGNUM
              // to avoid overflow when dividing by A[j][j].

              REC = (TJJ * BIGNUM) / XJ;
              if (CNORM[J] > ONE) {
                // Scale by 1/CNORM[j] to avoid overflow when
                // multiplying x[j] times column j.

                REC /= CNORM[J];
              }
              dscal(N, REC, X, 1);
              SCALE.value *= REC;
              XMAX *= REC;
            }
            X[J] /= TJJS;
            XJ = X[J].abs();
          } else {
            // A[j][j] = 0:  Set x[1:n] = 0, x[j] = 1, and
            // scale = 0, and compute a solution to A*x = 0.

            for (var I = 1; I <= N; I++) {
              X[I] = ZERO;
            }
            X[J] = ONE;
            XJ = ONE;
            SCALE.value = ZERO;
            XMAX = ZERO;
          }
        }

        // Scale x if necessary to avoid overflow when adding a
        // multiple of column j of A.

        if (XJ > ONE) {
          REC = ONE / XJ;
          if (CNORM[J] > (BIGNUM - XMAX) * REC) {
            // Scale x by 1/(2*abs(x[j])).

            REC *= HALF;
            dscal(N, REC, X, 1);
            SCALE.value *= REC;
          }
        } else if (XJ * CNORM[J] > (BIGNUM - XMAX)) {
          // Scale x by 1/2.

          dscal(N, HALF, X, 1);
          SCALE.value *= HALF;
        }

        if (UPPER) {
          if (J > 1) {
            // Compute the update
            // x[1:j-1] := x[1:j-1] - x[j] * A[1:j-1][j]

            daxpy(J - 1, -X[J] * TSCAL, A(1, J).asArray(), 1, X, 1);
            I = idamax(J - 1, X, 1);
            XMAX = X[I].abs();
          }
        } else {
          if (J < N) {
            // Compute the update
            // x[j+1:n] := x[j+1:n] - x[j] * A[j+1:n][j]

            daxpy(N - J, -X[J] * TSCAL, A(J + 1, J).asArray(), 1, X(J + 1), 1);
            I = J + idamax(N - J, X(J + 1), 1);
            XMAX = X[I].abs();
          }
        }
      }
    } else {
      // Solve A**T * x = b

      for (var J = JFIRST; JINC < 0 ? J >= JLAST : J <= JLAST; J += JINC) {
        // Compute x[j] = b(j) - sum A[k][j]*x[k].
        // k<>j

        XJ = X[J].abs();
        USCAL = TSCAL;
        REC = ONE / max(XMAX, ONE);
        if (CNORM[J] > (BIGNUM - XJ) * REC) {
          // If x[j] could overflow, scale x by 1/(2*XMAX).

          REC *= HALF;
          if (NOUNIT) {
            TJJS = A[J][J] * TSCAL;
          } else {
            TJJS = TSCAL;
          }
          TJJ = TJJS.abs();
          if (TJJ > ONE) {
            // Divide by A[j][j] when scaling x if A[j][j] > 1.

            REC = min(ONE, REC * TJJ);
            USCAL /= TJJS;
          }
          if (REC < ONE) {
            dscal(N, REC, X, 1);
            SCALE.value *= REC;
            XMAX *= REC;
          }
        }

        SUMJ = ZERO;
        if (USCAL == ONE) {
          // If the scaling needed for A in the dot product is 1,
          // call DDOT to perform the dot product.

          if (UPPER) {
            SUMJ = ddot(J - 1, A(1, J).asArray(), 1, X, 1);
          } else if (J < N) {
            SUMJ = ddot(N - J, A(J + 1, J).asArray(), 1, X(J + 1), 1);
          }
        } else {
          // Otherwise, use in-line code for the dot product.

          if (UPPER) {
            for (var I = 1; I <= J - 1; I++) {
              SUMJ += (A[I][J] * USCAL) * X[I];
            }
          } else if (J < N) {
            for (var I = J + 1; I <= N; I++) {
              SUMJ += (A[I][J] * USCAL) * X[I];
            }
          }
        }

        if (USCAL == TSCAL) {
          // Compute x[j] := ( x[j] - sumj ) / A[j][j] if 1/A[j][j]
          // was not used to scale the dotproduct.

          X[J] -= SUMJ;
          XJ = X[J].abs();
          var scale = true;
          if (NOUNIT) {
            TJJS = A[J][J] * TSCAL;
          } else {
            TJJS = TSCAL;
            if (TSCAL == ONE) {
              scale = false;
            }
          }

          if (scale) {
            // Compute x[j] /= A[j][j], scaling if necessary.

            TJJ = TJJS.abs();
            if (TJJ > SMLNUM) {
              // abs(A[j][j]) > SMLNUM:

              if (TJJ < ONE) {
                if (XJ > TJJ * BIGNUM) {
                  // Scale X by 1/abs(x[j]).

                  REC = ONE / XJ;
                  dscal(N, REC, X, 1);
                  SCALE.value *= REC;
                  XMAX *= REC;
                }
              }
              X[J] /= TJJS;
            } else if (TJJ > ZERO) {
              // 0 < abs(A[j][j]) <= SMLNUM:

              if (XJ > TJJ * BIGNUM) {
                // Scale x by (1/abs(x[j]))*abs(A[j][j])*BIGNUM.

                REC = (TJJ * BIGNUM) / XJ;
                dscal(N, REC, X, 1);
                SCALE.value *= REC;
                XMAX *= REC;
              }
              X[J] /= TJJS;
            } else {
              // A[j][j] = 0:  Set x[1:n] = 0, x[j] = 1, and
              // scale = 0, and compute a solution to A**T*x = 0.

              for (var I = 1; I <= N; I++) {
                X[I] = ZERO;
              }
              X[J] = ONE;
              SCALE.value = ZERO;
              XMAX = ZERO;
            }
          }
        } else {
          // Compute x[j] := x[j] / A[j][j]  - sumj if the dot
          // product has already been divided by 1/A[j][j].

          X[J] = X[J] / TJJS - SUMJ;
        }
        XMAX = max(XMAX, X[J].abs());
      }
    }
    SCALE.value /= TSCAL;
  }

  // Scale the column norms by 1/TSCAL for return.

  if (TSCAL != ONE) {
    dscal(N, ONE / TSCAL, CNORM, 1);
  }
}
