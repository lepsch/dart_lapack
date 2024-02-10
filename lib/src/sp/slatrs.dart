      void slatrs(UPLO, TRANS, DIAG, NORMIN, N, final Matrix<double> A, final int LDA, X, SCALE, CNORM, final Box<int> INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, NORMIN, TRANS, UPLO;
      int                INFO, LDA, N;
      double               SCALE;
      double               A( LDA, * ), CNORM( * ), X( * );
      // ..

      double               ZERO, HALF, ONE;
      const              ZERO = 0.0, HALF = 0.5, ONE = 1.0 ;
      bool               NOTRAN, NOUNIT, UPPER;
      int                I, IMAX, J, JFIRST, JINC, JLAST;
      double               BIGNUM, GROW, REC, SMLNUM, SUMJ, TJJ, TJJS, TMAX, TSCAL, USCAL, XBND, XJ, XMAX;
      double               WORK(1);
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ISAMAX;
      //- REAL               SASUM, SDOT, SLAMCH, SLANGE;
      // EXTERNAL lsame, ISAMAX, SASUM, SDOT, SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SSCAL, STRSV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      NOTRAN = lsame( TRANS, 'N' );
      NOUNIT = lsame( DIAG, 'N' );

      // Test the input parameters.

      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( !NOTRAN && !lsame( TRANS, 'T' ) && !lsame( TRANS, 'C' ) ) {
         INFO = -2;
      } else if ( !NOUNIT && !lsame( DIAG, 'U' ) ) {
         INFO = -3;
      } else if ( !lsame( NORMIN, 'Y' ) && !lsame( NORMIN, 'N' ) ) {
         INFO = -4;
      } else if ( N < 0 ) {
         INFO = -5;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -7;
      }
      if ( INFO != 0 ) {
         xerbla('SLATRS', -INFO );
         return;
      }

      // Quick return if possible

      SCALE = ONE;
      if (N == 0) return;

      // Determine machine dependent parameters to control overflow.

      SMLNUM = SLAMCH( 'Safe minimum' ) / SLAMCH( 'Precision' );
      BIGNUM = ONE / SMLNUM;

      if ( lsame( NORMIN, 'N' ) ) {

         // Compute the 1-norm of each column, not including the diagonal.

         if ( UPPER ) {

            // A is upper triangular.

            for (J = 1; J <= N; J++) { // 10
               CNORM[J] = SASUM( J-1, A( 1, J ), 1 );
            } // 10
         } else {

            // A is lower triangular.

            for (J = 1; J <= N - 1; J++) { // 20
               CNORM[J] = SASUM( N-J, A( J+1, J ), 1 );
            } // 20
            CNORM[N] = ZERO;
         }
      }

      // Scale the column norms by TSCAL if the maximum element in CNORM is
      // greater than BIGNUM.

      IMAX = ISAMAX( N, CNORM, 1 );
      TMAX = CNORM( IMAX );
      if ( TMAX <= BIGNUM ) {
         TSCAL = ONE;
      } else {

         // Avoid NaN generation if entries in CNORM exceed the
         // overflow threshold

         if ( TMAX <= SLAMCH('Overflow') ) {
            // Case 1: All entries in CNORM are valid floating-point numbers
            TSCAL = ONE / ( SMLNUM*TMAX );
            sscal(N, TSCAL, CNORM, 1 );
         } else {
            // Case 2: At least one column norm of A cannot be represented
            // as floating-point number. Find the offdiagonal entry A( I, J )
            // with the largest absolute value. If this entry is not +/- Infinity,
            // use this value as TSCAL.
            TMAX = ZERO;
            if ( UPPER ) {

               // A is upper triangular.

               for (J = 2; J <= N; J++) {
                  TMAX = max( SLANGE( 'M', J-1, 1, A( 1, J ), 1, WORK ), TMAX );
               }
            } else {

               // A is lower triangular.

               for (J = 1; J <= N - 1; J++) {
                  TMAX = max( SLANGE( 'M', N-J, 1, A( J+1, J ), 1, WORK ), TMAX );
               }
            }

            if ( TMAX <= SLAMCH('Overflow') ) {
               TSCAL = ONE / ( SMLNUM*TMAX );
               for (J = 1; J <= N; J++) {
                  if ( CNORM( J ) <= SLAMCH('Overflow') ) {
                     CNORM[J] = CNORM( J )*TSCAL;
                  } else {
                     // Recompute the 1-norm without introducing Infinity
                     // in the summation
                     CNORM[J] = ZERO;
                     if ( UPPER ) {
                        for (I = 1; I <= J - 1; I++) {
                           CNORM[J] = CNORM( J ) + TSCAL * ( A( I, J ) ).abs();
                        }
                     } else {
                        for (I = J + 1; I <= N; I++) {
                           CNORM[J] = CNORM( J ) + TSCAL * ( A( I, J ) ).abs();
                        }
                     }
                  }
               }
            } else {
               // At least one entry of A is not a valid floating-point entry.
               // Rely on TRSV to propagate Inf and NaN.
               strsv(UPLO, TRANS, DIAG, N, A, LDA, X, 1 );
               return;
            }
         }
      }

      // Compute a bound on the computed solution vector to see if the
      // Level 2 BLAS routine STRSV can be used.

      J = ISAMAX( N, X, 1 );
      XMAX = ( X( J ) ).abs();
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
            for (J = JFIRST; JINC < 0 ? J >= JLAST : J <= JLAST; J += JINC) { // 30

               // Exit the loop if the growth factor is too small.

               if (GROW <= SMLNUM) GO TO 50;

               // M(j) = G(j-1) / abs(A(j,j))

               TJJ = ( A( J, J ) ).abs();
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
            for (J = JFIRST; JINC < 0 ? J >= JLAST : J <= JLAST; J += JINC) { // 40

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
            for (J = JFIRST; JINC < 0 ? J >= JLAST : J <= JLAST; J += JINC) { // 60

               // Exit the loop if the growth factor is too small.

               if (GROW <= SMLNUM) GO TO 80;

               // G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )

               XJ = ONE + CNORM( J );
               GROW = min( GROW, XBND / XJ );

               // M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))

               TJJ = ( A( J, J ) ).abs();
               if (XJ > TJJ) XBND = XBND*( TJJ / XJ );
            } // 60
            GROW = min( GROW, XBND );
         } else {

            // A is unit triangular.

            // Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.

            GROW = min( ONE, ONE / max( XBND, SMLNUM ) );
            for (J = JFIRST; JINC < 0 ? J >= JLAST : J <= JLAST; J += JINC) { // 70

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

         strsv(UPLO, TRANS, DIAG, N, A, LDA, X, 1 );
      } else {

         // Use a Level 1 BLAS solve, scaling intermediate results.

         if ( XMAX > BIGNUM ) {

            // Scale X so that its components are less than or equal to
            // BIGNUM in absolute value.

            SCALE = BIGNUM / XMAX;
            sscal(N, SCALE, X, 1 );
            XMAX = BIGNUM;
         }

         if ( NOTRAN ) {

            // Solve A * x = b

            for (J = JFIRST; JINC < 0 ? J >= JLAST : J <= JLAST; J += JINC) { // 100

               // Compute x(j) = b(j) / A(j,j), scaling x if necessary.

               XJ = ( X( J ) ).abs();
               if ( NOUNIT ) {
                  TJJS = A( J, J )*TSCAL;
               } else {
                  TJJS = TSCAL;
                  if (TSCAL == ONE) GO TO 95;
               }
                  TJJ = ( TJJS ).abs();
                  if ( TJJ > SMLNUM ) {

                     // abs(A(j,j)) > SMLNUM:

                     if ( TJJ < ONE ) {
                        if ( XJ > TJJ*BIGNUM ) {

                           // Scale x by 1/b(j).

                           REC = ONE / XJ;
                           sscal(N, REC, X, 1 );
                           SCALE = SCALE*REC;
                           XMAX = XMAX*REC;
                        }
                     }
                     X[J] = X( J ) / TJJS;
                     XJ = ( X( J ) ).abs();
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
                        sscal(N, REC, X, 1 );
                        SCALE = SCALE*REC;
                        XMAX = XMAX*REC;
                     }
                     X[J] = X( J ) / TJJS;
                     XJ = ( X( J ) ).abs();
                  } else {

                     // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                     // scale = 0, and compute a solution to A*x = 0.

                     for (I = 1; I <= N; I++) { // 90
                        X[I] = ZERO;
                     } // 90
                     X[J] = ONE;
                     XJ = ONE;
                     SCALE = ZERO;
                     XMAX = ZERO;
                  }
               } // 95

               // Scale x if necessary to avoid overflow when adding a
               // multiple of column j of A.

               if ( XJ > ONE ) {
                  REC = ONE / XJ;
                  if ( CNORM( J ) > ( BIGNUM-XMAX )*REC ) {

                     // Scale x by 1/(2*abs(x(j))).

                     REC = REC*HALF;
                     sscal(N, REC, X, 1 );
                     SCALE = SCALE*REC;
                  }
               } else if ( XJ*CNORM( J ) > ( BIGNUM-XMAX ) ) {

                  // Scale x by 1/2.

                  sscal(N, HALF, X, 1 );
                  SCALE = SCALE*HALF;
               }

               if ( UPPER ) {
                  if ( J > 1 ) {

                     // Compute the update
                     //    x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)

                     saxpy(J-1, -X( J )*TSCAL, A( 1, J ), 1, X, 1 );
                     I = ISAMAX( J-1, X, 1 );
                     XMAX = ( X( I ) ).abs();
                  }
               } else {
                  if ( J < N ) {

                     // Compute the update
                     //    x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)

                     saxpy(N-J, -X( J )*TSCAL, A( J+1, J ), 1, X( J+1 ), 1 );
                     I = J + ISAMAX( N-J, X( J+1 ), 1 );
                     XMAX = ( X( I ) ).abs();
                  }
               }
            } // 100

         } else {

            // Solve A**T * x = b

            for (J = JFIRST; JINC < 0 ? J >= JLAST : J <= JLAST; J += JINC) { // 140

               // Compute x(j) = b(j) - sum A(k,j)*x(k).
               //                       k<>j

               XJ = ( X( J ) ).abs();
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
                     TJJ = ( TJJS ).abs();
                     if ( TJJ > ONE ) {

                        // Divide by A(j,j) when scaling x if A(j,j) > 1.

                        REC = min( ONE, REC*TJJ );
                        USCAL = USCAL / TJJS;
                     }
                  if ( REC < ONE ) {
                     sscal(N, REC, X, 1 );
                     SCALE = SCALE*REC;
                     XMAX = XMAX*REC;
                  }
               }

               SUMJ = ZERO;
               if ( USCAL == ONE ) {

                  // If the scaling needed for A in the dot product is 1,
                  // call SDOT to perform the dot product.

                  if ( UPPER ) {
                     SUMJ = SDOT( J-1, A( 1, J ), 1, X, 1 );
                  } else if ( J < N ) {
                     SUMJ = SDOT( N-J, A( J+1, J ), 1, X( J+1 ), 1 );
                  }
               } else {

                  // Otherwise, use in-line code for the dot product.

                  if ( UPPER ) {
                     for (I = 1; I <= J - 1; I++) { // 110
                        SUMJ = SUMJ + ( A( I, J )*USCAL )*X( I );
                     } // 110
                  } else if ( J < N ) {
                     for (I = J + 1; I <= N; I++) { // 120
                        SUMJ = SUMJ + ( A( I, J )*USCAL )*X( I );
                     } // 120
                  }
               }

               if ( USCAL == TSCAL ) {

                  // Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j)
                  // was not used to scale the dotproduct.

                  X[J] = X( J ) - SUMJ;
                  XJ = ( X( J ) ).abs();
                  if ( NOUNIT ) {
                     TJJS = A( J, J )*TSCAL;
                  } else {
                     TJJS = TSCAL;
                     if (TSCAL == ONE) GO TO 135;
                  }

                     // Compute x(j) = x(j) / A(j,j), scaling if necessary.

                     TJJ = ( TJJS ).abs();
                     if ( TJJ > SMLNUM ) {

                        // abs(A(j,j)) > SMLNUM:

                        if ( TJJ < ONE ) {
                           if ( XJ > TJJ*BIGNUM ) {

                              // Scale X by 1/abs(x(j)).

                              REC = ONE / XJ;
                              sscal(N, REC, X, 1 );
                              SCALE = SCALE*REC;
                              XMAX = XMAX*REC;
                           }
                        }
                        X[J] = X( J ) / TJJS;
                     } else if ( TJJ > ZERO ) {

                        // 0 < abs(A(j,j)) <= SMLNUM:

                        if ( XJ > TJJ*BIGNUM ) {

                           // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.

                           REC = ( TJJ*BIGNUM ) / XJ;
                           sscal(N, REC, X, 1 );
                           SCALE = SCALE*REC;
                           XMAX = XMAX*REC;
                        }
                        X[J] = X( J ) / TJJS;
                     } else {

                        // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                        // scale = 0, and compute a solution to A**T*x = 0.

                        for (I = 1; I <= N; I++) { // 130
                           X[I] = ZERO;
                        } // 130
                        X[J] = ONE;
                        SCALE = ZERO;
                        XMAX = ZERO;
                     }
                  } // 135
               } else {

                  // Compute x(j) := x(j) / A(j,j)  - sumj if the dot
                  // product has already been divided by 1/A(j,j).

                  X[J] = X( J ) / TJJS - SUMJ;
               }
               XMAX = max( XMAX, ( X( J ) ).abs() );
            } // 140
         }
         SCALE = SCALE / TSCAL;
      }

      // Scale the column norms by 1/TSCAL for return.

      if ( TSCAL != ONE ) {
         sscal(N, ONE / TSCAL, CNORM, 1 );
      }

      }
