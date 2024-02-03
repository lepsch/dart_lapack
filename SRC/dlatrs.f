      SUBROUTINE DLATRS( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE, CNORM, INFO );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, NORMIN, TRANS, UPLO;
      int                INFO, LDA, N;
      double             SCALE;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), CNORM( * ), X( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, HALF, ONE;
      const              ZERO = 0.0, HALF = 0.5, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRAN, NOUNIT, UPPER;
      int                I, IMAX, J, JFIRST, JINC, JLAST;
      double             BIGNUM, GROW, REC, SMLNUM, SUMJ, TJJ, TJJS, TMAX, TSCAL, USCAL, XBND, XJ, XMAX;
      // ..
      // .. Local Arrays ..
      double             WORK(1);
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IDAMAX;
      double             DASUM, DDOT, DLAMCH, DLANGE;
      // EXTERNAL LSAME, IDAMAX, DASUM, DDOT, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DSCAL, DTRSV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      NOTRAN = LSAME( TRANS, 'N' );
      NOUNIT = LSAME( DIAG, 'N' );

      // Test the input parameters.

      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( !NOTRAN && !LSAME( TRANS, 'T' ) && !LSAME( TRANS, 'C' ) ) {
         INFO = -2;
      } else if ( !NOUNIT && !LSAME( DIAG, 'U' ) ) {
         INFO = -3;
      } else if ( !LSAME( NORMIN, 'Y' ) && !LSAME( NORMIN, 'N' ) ) {
         INFO = -4;
      } else if ( N < 0 ) {
         INFO = -5;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -7;
      }
      if ( INFO != 0 ) {
         xerbla('DLATRS', -INFO );
         return;
      }

      // Quick return if possible

      SCALE = ONE;
      if (N == 0) return;

      // Determine machine dependent parameters to control overflow.

      SMLNUM = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' );
      BIGNUM = ONE / SMLNUM;

      if ( LSAME( NORMIN, 'N' ) ) {

         // Compute the 1-norm of each column, not including the diagonal.

         if ( UPPER ) {

            // A is upper triangular.

            for (J = 1; J <= N; J++) { // 10
               CNORM( J ) = DASUM( J-1, A( 1, J ), 1 );
            } // 10
         } else {

            // A is lower triangular.

            for (J = 1; J <= N - 1; J++) { // 20
               CNORM( J ) = DASUM( N-J, A( J+1, J ), 1 );
            } // 20
            CNORM( N ) = ZERO;
         }
      }

      // Scale the column norms by TSCAL if the maximum element in CNORM is
      // greater than BIGNUM.

      IMAX = IDAMAX( N, CNORM, 1 );
      TMAX = CNORM( IMAX );
      if ( TMAX <= BIGNUM ) {
         TSCAL = ONE;
      } else {

         // Avoid NaN generation if entries in CNORM exceed the
         // overflow threshold

         if ( TMAX <= DLAMCH('Overflow') ) {
            // Case 1: All entries in CNORM are valid floating-point numbers
            TSCAL = ONE / ( SMLNUM*TMAX );
            dscal(N, TSCAL, CNORM, 1 );
         } else {
            // Case 2: At least one column norm of A cannot be represented
            // as floating-point number. Find the offdiagonal entry A( I, J )
            // with the largest absolute value. If this entry is not +/- Infinity,
            // use this value as TSCAL.
            TMAX = ZERO;
            if ( UPPER ) {

               // A is upper triangular.

               for (J = 2; J <= N; J++) {
                  TMAX = max( DLANGE( 'M', J-1, 1, A( 1, J ), 1, WORK ), TMAX );
               }
            } else {

               // A is lower triangular.

               for (J = 1; J <= N - 1; J++) {
                  TMAX = max( DLANGE( 'M', N-J, 1, A( J+1, J ), 1, WORK ), TMAX );
               }
            }

            if ( TMAX <= DLAMCH('Overflow') ) {
               TSCAL = ONE / ( SMLNUM*TMAX );
               for (J = 1; J <= N; J++) {
                  if ( CNORM( J ) <= DLAMCH('Overflow') ) {
                     CNORM( J ) = CNORM( J )*TSCAL;
                  } else {
                     // Recompute the 1-norm without introducing Infinity
                     // in the summation
                     CNORM( J ) = ZERO;
                     if ( UPPER ) {
                        for (I = 1; I <= J - 1; I++) {
                           CNORM( J ) = CNORM( J ) + TSCAL * ABS( A( I, J ) );
                        }
                     } else {
                        for (I = J + 1; I <= N; I++) {
                           CNORM( J ) = CNORM( J ) + TSCAL * ABS( A( I, J ) );
                        }
                     }
                  }
               }
            } else {
               // At least one entry of A is not a valid floating-point entry.
               // Rely on TRSV to propagate Inf and NaN.
               dtrsv(UPLO, TRANS, DIAG, N, A, LDA, X, 1 );
               return;
            }
         }
      }

      // Compute a bound on the computed solution vector to see if the
      // Level 2 BLAS routine DTRSV can be used.

      J = IDAMAX( N, X, 1 );
      XMAX = ABS( X( J ) );
      XBND = XMAX;
      if ( NOTRAN ) {

         // Compute the growth in A * x = b.

         if ( UPPER ) {
            JFIRST = N;
            JLAST = 1;
            JINC = -1;
         } else {
            JFIRST = 1;
            JLAST = N;
            JINC = 1;
         }

         if ( TSCAL != ONE ) {
            GROW = ZERO;
            GO TO 50;
         }

         if ( NOUNIT ) {

            // A is non-unit triangular.

            // Compute GROW = 1/G(j) and XBND = 1/M(j).
            // Initially, G(0) = max{x(i), i=1,...,n}.

            GROW = ONE / max( XBND, SMLNUM );
            XBND = GROW;
            DO 30 J = JFIRST, JLAST, JINC;

               // Exit the loop if the growth factor is too small.

               if (GROW <= SMLNUM) GO TO 50;

               // M(j) = G(j-1) / abs(A(j,j))

               TJJ = ABS( A( J, J ) );
               XBND = min( XBND, min( ONE, TJJ )*GROW );
               if ( TJJ+CNORM( J ) >= SMLNUM ) {

                  // G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )

                  GROW = GROW*( TJJ / ( TJJ+CNORM( J ) ) );
               } else {

                  // G(j) could overflow, set GROW to 0.

                  GROW = ZERO;
               }
            } // 30
            GROW = XBND;
         } else {

            // A is unit triangular.

            // Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.

            GROW = min( ONE, ONE / max( XBND, SMLNUM ) );
            DO 40 J = JFIRST, JLAST, JINC;

               // Exit the loop if the growth factor is too small.

               if (GROW <= SMLNUM) GO TO 50;

               // G(j) = G(j-1)*( 1 + CNORM(j) )

               GROW = GROW*( ONE / ( ONE+CNORM( J ) ) );
            } // 40
         }
         } // 50

      } else {

         // Compute the growth in A**T * x = b.

         if ( UPPER ) {
            JFIRST = 1;
            JLAST = N;
            JINC = 1;
         } else {
            JFIRST = N;
            JLAST = 1;
            JINC = -1;
         }

         if ( TSCAL != ONE ) {
            GROW = ZERO;
            GO TO 80;
         }

         if ( NOUNIT ) {

            // A is non-unit triangular.

            // Compute GROW = 1/G(j) and XBND = 1/M(j).
            // Initially, M(0) = max{x(i), i=1,...,n}.

            GROW = ONE / max( XBND, SMLNUM );
            XBND = GROW;
            DO 60 J = JFIRST, JLAST, JINC;

               // Exit the loop if the growth factor is too small.

               if (GROW <= SMLNUM) GO TO 80;

               // G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )

               XJ = ONE + CNORM( J );
               GROW = min( GROW, XBND / XJ );

               // M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))

               TJJ = ABS( A( J, J ) );
               if (XJ > TJJ) XBND = XBND*( TJJ / XJ );
            } // 60
            GROW = min( GROW, XBND );
         } else {

            // A is unit triangular.

            // Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.

            GROW = min( ONE, ONE / max( XBND, SMLNUM ) );
            DO 70 J = JFIRST, JLAST, JINC;

               // Exit the loop if the growth factor is too small.

               if (GROW <= SMLNUM) GO TO 80;

               // G(j) = ( 1 + CNORM(j) )*G(j-1)

               XJ = ONE + CNORM( J );
               GROW = GROW / XJ;
            } // 70
         }
         } // 80
      }

      if ( ( GROW*TSCAL ) > SMLNUM ) {

         // Use the Level 2 BLAS solve if the reciprocal of the bound on
         // elements of X is not too small.

         dtrsv(UPLO, TRANS, DIAG, N, A, LDA, X, 1 );
      } else {

         // Use a Level 1 BLAS solve, scaling intermediate results.

         if ( XMAX > BIGNUM ) {

            // Scale X so that its components are less than or equal to
            // BIGNUM in absolute value.

            SCALE = BIGNUM / XMAX;
            dscal(N, SCALE, X, 1 );
            XMAX = BIGNUM;
         }

         if ( NOTRAN ) {

            // Solve A * x = b

            DO 110 J = JFIRST, JLAST, JINC;

               // Compute x(j) = b(j) / A(j,j), scaling x if necessary.

               XJ = ABS( X( J ) );
               if ( NOUNIT ) {
                  TJJS = A( J, J )*TSCAL;
               } else {
                  TJJS = TSCAL;
                  if (TSCAL == ONE) GO TO 100;
               }
               TJJ = ABS( TJJS );
               if ( TJJ > SMLNUM ) {

                     // abs(A(j,j)) > SMLNUM:

                  if ( TJJ < ONE ) {
                     if ( XJ > TJJ*BIGNUM ) {

                           // Scale x by 1/b(j).

                        REC = ONE / XJ;
                        dscal(N, REC, X, 1 );
                        SCALE = SCALE*REC;
                        XMAX = XMAX*REC;
                     }
                  }
                  X( J ) = X( J ) / TJJS;
                  XJ = ABS( X( J ) );
               } else if ( TJJ > ZERO ) {

                     // 0 < abs(A(j,j)) <= SMLNUM:

                  if ( XJ > TJJ*BIGNUM ) {

                        // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
                        // to avoid overflow when dividing by A(j,j).

                     REC = ( TJJ*BIGNUM ) / XJ;
                     if ( CNORM( J ) > ONE ) {

                           // Scale by 1/CNORM(j) to avoid overflow when
                           // multiplying x(j) times column j.

                        REC = REC / CNORM( J );
                     }
                     dscal(N, REC, X, 1 );
                     SCALE = SCALE*REC;
                     XMAX = XMAX*REC;
                  }
                  X( J ) = X( J ) / TJJS;
                  XJ = ABS( X( J ) );
               } else {

                     // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                     // scale = 0, and compute a solution to A*x = 0.

                  for (I = 1; I <= N; I++) { // 90
                     X( I ) = ZERO;
                  } // 90
                  X( J ) = ONE;
                  XJ = ONE;
                  SCALE = ZERO;
                  XMAX = ZERO;
               }
               } // 100

               // Scale x if necessary to avoid overflow when adding a
               // multiple of column j of A.

               if ( XJ > ONE ) {
                  REC = ONE / XJ;
                  if ( CNORM( J ) > ( BIGNUM-XMAX )*REC ) {

                     // Scale x by 1/(2*abs(x(j))).

                     REC = REC*HALF;
                     dscal(N, REC, X, 1 );
                     SCALE = SCALE*REC;
                  }
               } else if ( XJ*CNORM( J ) > ( BIGNUM-XMAX ) ) {

                  // Scale x by 1/2.

                  dscal(N, HALF, X, 1 );
                  SCALE = SCALE*HALF;
               }

               if ( UPPER ) {
                  if ( J > 1 ) {

                     // Compute the update
                        // x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)

                     daxpy(J-1, -X( J )*TSCAL, A( 1, J ), 1, X, 1 );
                     I = IDAMAX( J-1, X, 1 );
                     XMAX = ABS( X( I ) );
                  }
               } else {
                  if ( J < N ) {

                     // Compute the update
                        // x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)

                     daxpy(N-J, -X( J )*TSCAL, A( J+1, J ), 1, X( J+1 ), 1 );
                     I = J + IDAMAX( N-J, X( J+1 ), 1 );
                     XMAX = ABS( X( I ) );
                  }
               }
            } // 110

         } else {

            // Solve A**T * x = b

            DO 160 J = JFIRST, JLAST, JINC;

               // Compute x(j) = b(j) - sum A(k,j)*x(k).
                                     // k<>j

               XJ = ABS( X( J ) );
               USCAL = TSCAL;
               REC = ONE / max( XMAX, ONE );
               if ( CNORM( J ) > ( BIGNUM-XJ )*REC ) {

                  // If x(j) could overflow, scale x by 1/(2*XMAX).

                  REC = REC*HALF;
                  if ( NOUNIT ) {
                     TJJS = A( J, J )*TSCAL;
                  } else {
                     TJJS = TSCAL;
                  }
                  TJJ = ABS( TJJS );
                  if ( TJJ > ONE ) {

                        // Divide by A(j,j) when scaling x if A(j,j) > 1.

                     REC = min( ONE, REC*TJJ );
                     USCAL = USCAL / TJJS;
                  }
                  if ( REC < ONE ) {
                     dscal(N, REC, X, 1 );
                     SCALE = SCALE*REC;
                     XMAX = XMAX*REC;
                  }
               }

               SUMJ = ZERO;
               if ( USCAL == ONE ) {

                  // If the scaling needed for A in the dot product is 1,
                  // call DDOT to perform the dot product.

                  if ( UPPER ) {
                     SUMJ = DDOT( J-1, A( 1, J ), 1, X, 1 );
                  } else if ( J < N ) {
                     SUMJ = DDOT( N-J, A( J+1, J ), 1, X( J+1 ), 1 );
                  }
               } else {

                  // Otherwise, use in-line code for the dot product.

                  if ( UPPER ) {
                     for (I = 1; I <= J - 1; I++) { // 120
                        SUMJ = SUMJ + ( A( I, J )*USCAL )*X( I );
                     } // 120
                  } else if ( J < N ) {
                     for (I = J + 1; I <= N; I++) { // 130
                        SUMJ = SUMJ + ( A( I, J )*USCAL )*X( I );
                     } // 130
                  }
               }

               if ( USCAL == TSCAL ) {

                  // Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j)
                  // was not used to scale the dotproduct.

                  X( J ) = X( J ) - SUMJ;
                  XJ = ABS( X( J ) );
                  if ( NOUNIT ) {
                     TJJS = A( J, J )*TSCAL;
                  } else {
                     TJJS = TSCAL;
                     if (TSCAL == ONE) GO TO 150;
                  }

                     // Compute x(j) = x(j) / A(j,j), scaling if necessary.

                  TJJ = ABS( TJJS );
                  if ( TJJ > SMLNUM ) {

                        // abs(A(j,j)) > SMLNUM:

                     if ( TJJ < ONE ) {
                        if ( XJ > TJJ*BIGNUM ) {

                              // Scale X by 1/abs(x(j)).

                           REC = ONE / XJ;
                           dscal(N, REC, X, 1 );
                           SCALE = SCALE*REC;
                           XMAX = XMAX*REC;
                        }
                     }
                     X( J ) = X( J ) / TJJS;
                  } else if ( TJJ > ZERO ) {

                        // 0 < abs(A(j,j)) <= SMLNUM:

                     if ( XJ > TJJ*BIGNUM ) {

                           // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.

                        REC = ( TJJ*BIGNUM ) / XJ;
                        dscal(N, REC, X, 1 );
                        SCALE = SCALE*REC;
                        XMAX = XMAX*REC;
                     }
                     X( J ) = X( J ) / TJJS;
                  } else {

                        // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                        // scale = 0, and compute a solution to A**T*x = 0.

                     for (I = 1; I <= N; I++) { // 140
                        X( I ) = ZERO;
                     } // 140
                     X( J ) = ONE;
                     SCALE = ZERO;
                     XMAX = ZERO;
                  }
                  } // 150
               } else {

                  // Compute x(j) := x(j) / A(j,j)  - sumj if the dot
                  // product has already been divided by 1/A(j,j).

                  X( J ) = X( J ) / TJJS - SUMJ;
               }
               XMAX = max( XMAX, ABS( X( J ) ) );
            } // 160
         }
         SCALE = SCALE / TSCAL;
      }

      // Scale the column norms by 1/TSCAL for return.

      if ( TSCAL != ONE ) {
         dscal(N, ONE / TSCAL, CNORM, 1 );
      }

      return;

      // End of DLATRS

      }
