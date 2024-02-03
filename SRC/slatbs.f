      SUBROUTINE SLATBS( UPLO, TRANS, DIAG, NORMIN, N, KD, AB, LDAB, X, SCALE, CNORM, INFO );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, NORMIN, TRANS, UPLO;
      int                INFO, KD, LDAB, N;
      REAL               SCALE;
      // ..
      // .. Array Arguments ..
      REAL               AB( LDAB, * ), CNORM( * ), X( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, HALF, ONE;
      const              ZERO = 0.0, HALF = 0.5, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRAN, NOUNIT, UPPER;
      int                I, IMAX, J, JFIRST, JINC, JLAST, JLEN, MAIND;
      REAL               BIGNUM, GROW, REC, SMLNUM, SUMJ, TJJ, TJJS, TMAX, TSCAL, USCAL, XBND, XJ, XMAX;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ISAMAX;
      REAL               SASUM, SDOT, SLAMCH;
      // EXTERNAL LSAME, ISAMAX, SASUM, SDOT, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SSCAL, STBSV, XERBLA
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
      } else if ( KD < 0 ) {
         INFO = -6;
      } else if ( LDAB < KD+1 ) {
         INFO = -8;
      }
      if ( INFO != 0 ) {
         xerbla('SLATBS', -INFO );
         return;
      }

      // Quick return if possible

      SCALE = ONE;
      if (N == 0) RETURN;

      // Determine machine dependent parameters to control overflow.

      SMLNUM = SLAMCH( 'Safe minimum' ) / SLAMCH( 'Precision' );
      BIGNUM = ONE / SMLNUM;

      if ( LSAME( NORMIN, 'N' ) ) {

         // Compute the 1-norm of each column, not including the diagonal.

         if ( UPPER ) {

            // A is upper triangular.

            for (J = 1; J <= N; J++) { // 10
               JLEN = MIN( KD, J-1 );
               CNORM( J ) = SASUM( JLEN, AB( KD+1-JLEN, J ), 1 );
            } // 10
         } else {

            // A is lower triangular.

            for (J = 1; J <= N; J++) { // 20
               JLEN = MIN( KD, N-J );
               if ( JLEN > 0 ) {
                  CNORM( J ) = SASUM( JLEN, AB( 2, J ), 1 );
               } else {
                  CNORM( J ) = ZERO;
               }
            } // 20
         }
      }

      // Scale the column norms by TSCAL if the maximum element in CNORM is
      // greater than BIGNUM.

      IMAX = ISAMAX( N, CNORM, 1 );
      TMAX = CNORM( IMAX );
      if ( TMAX <= BIGNUM ) {
         TSCAL = ONE;
      } else {
         TSCAL = ONE / ( SMLNUM*TMAX );
         sscal(N, TSCAL, CNORM, 1 );
      }

      // Compute a bound on the computed solution vector to see if the
      // Level 2 BLAS routine STBSV can be used.

      J = ISAMAX( N, X, 1 );
      XMAX = ABS( X( J ) );
      XBND = XMAX;
      if ( NOTRAN ) {

         // Compute the growth in A * x = b.

         if ( UPPER ) {
            JFIRST = N;
            JLAST = 1;
            JINC = -1;
            MAIND = KD + 1;
         } else {
            JFIRST = 1;
            JLAST = N;
            JINC = 1;
            MAIND = 1;
         }

         if ( TSCAL != ONE ) {
            GROW = ZERO;
            GO TO 50;
         }

         if ( NOUNIT ) {

            // A is non-unit triangular.

            // Compute GROW = 1/G(j) and XBND = 1/M(j).
            // Initially, G(0) = max{x(i), i=1,...,n}.

            GROW = ONE / MAX( XBND, SMLNUM );
            XBND = GROW;
            DO 30 J = JFIRST, JLAST, JINC;

               // Exit the loop if the growth factor is too small.

               if (GROW <= SMLNUM) GO TO 50;

               // M(j) = G(j-1) / abs(A(j,j))

               TJJ = ABS( AB( MAIND, J ) );
               XBND = MIN( XBND, MIN( ONE, TJJ )*GROW );
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

            GROW = MIN( ONE, ONE / MAX( XBND, SMLNUM ) );
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
            MAIND = KD + 1;
         } else {
            JFIRST = N;
            JLAST = 1;
            JINC = -1;
            MAIND = 1;
         }

         if ( TSCAL != ONE ) {
            GROW = ZERO;
            GO TO 80;
         }

         if ( NOUNIT ) {

            // A is non-unit triangular.

            // Compute GROW = 1/G(j) and XBND = 1/M(j).
            // Initially, M(0) = max{x(i), i=1,...,n}.

            GROW = ONE / MAX( XBND, SMLNUM );
            XBND = GROW;
            DO 60 J = JFIRST, JLAST, JINC;

               // Exit the loop if the growth factor is too small.

               if (GROW <= SMLNUM) GO TO 80;

               // G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )

               XJ = ONE + CNORM( J );
               GROW = MIN( GROW, XBND / XJ );

               // M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))

               TJJ = ABS( AB( MAIND, J ) );
               if (XJ > TJJ) XBND = XBND*( TJJ / XJ );
            } // 60
            GROW = MIN( GROW, XBND );
         } else {

            // A is unit triangular.

            // Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.

            GROW = MIN( ONE, ONE / MAX( XBND, SMLNUM ) );
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

         stbsv(UPLO, TRANS, DIAG, N, KD, AB, LDAB, X, 1 );
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

            DO 100 J = JFIRST, JLAST, JINC;

               // Compute x(j) = b(j) / A(j,j), scaling x if necessary.

               XJ = ABS( X( J ) );
               if ( NOUNIT ) {
                  TJJS = AB( MAIND, J )*TSCAL;
               } else {
                  TJJS = TSCAL;
                  if (TSCAL == ONE) GO TO 95;
               }
                  TJJ = ABS( TJJS );
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
                        sscal(N, REC, X, 1 );
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
                        // x(max(1,j-kd):j-1) := x(max(1,j-kd):j-1) -
                                              // x(j)* A(max(1,j-kd):j-1,j)

                     JLEN = MIN( KD, J-1 );
                     saxpy(JLEN, -X( J )*TSCAL, AB( KD+1-JLEN, J ), 1, X( J-JLEN ), 1 );
                     I = ISAMAX( J-1, X, 1 );
                     XMAX = ABS( X( I ) );
                  }
               } else if ( J < N ) {

                  // Compute the update
                     // x(j+1:min(j+kd,n)) := x(j+1:min(j+kd,n)) -
                                           // x(j) * A(j+1:min(j+kd,n),j)

                  JLEN = MIN( KD, N-J );
                  if (JLEN > 0) CALL SAXPY( JLEN, -X( J )*TSCAL, AB( 2, J ), 1, X( J+1 ), 1 );
                  I = J + ISAMAX( N-J, X( J+1 ), 1 );
                  XMAX = ABS( X( I ) );
               }
            } // 100

         } else {

            // Solve A**T * x = b

            DO 140 J = JFIRST, JLAST, JINC;

               // Compute x(j) = b(j) - sum A(k,j)*x(k).
                                     // k<>j

               XJ = ABS( X( J ) );
               USCAL = TSCAL;
               REC = ONE / MAX( XMAX, ONE );
               if ( CNORM( J ) > ( BIGNUM-XJ )*REC ) {

                  // If x(j) could overflow, scale x by 1/(2*XMAX).

                  REC = REC*HALF;
                  if ( NOUNIT ) {
                     TJJS = AB( MAIND, J )*TSCAL;
                  } else {
                     TJJS = TSCAL;
                  }
                     TJJ = ABS( TJJS );
                     if ( TJJ > ONE ) {

                        // Divide by A(j,j) when scaling x if A(j,j) > 1.

                        REC = MIN( ONE, REC*TJJ );
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
                     JLEN = MIN( KD, J-1 );
                     SUMJ = SDOT( JLEN, AB( KD+1-JLEN, J ), 1, X( J-JLEN ), 1 );
                  } else {
                     JLEN = MIN( KD, N-J );
                     if (JLEN > 0) SUMJ = SDOT( JLEN, AB( 2, J ), 1, X( J+1 ), 1 );
                  }
               } else {

                  // Otherwise, use in-line code for the dot product.

                  if ( UPPER ) {
                     JLEN = MIN( KD, J-1 );
                     for (I = 1; I <= JLEN; I++) { // 110
                        SUMJ = SUMJ + ( AB( KD+I-JLEN, J )*USCAL )* X( J-JLEN-1+I );
                     } // 110
                  } else {
                     JLEN = MIN( KD, N-J );
                     for (I = 1; I <= JLEN; I++) { // 120
                        SUMJ = SUMJ + ( AB( I+1, J )*USCAL )*X( J+I );
                     } // 120
                  }
               }

               if ( USCAL == TSCAL ) {

                  // Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j)
                  // was not used to scale the dotproduct.

                  X( J ) = X( J ) - SUMJ;
                  XJ = ABS( X( J ) );
                  if ( NOUNIT ) {

                     // Compute x(j) = x(j) / A(j,j), scaling if necessary.

                     TJJS = AB( MAIND, J )*TSCAL;
                  } else {
                     TJJS = TSCAL;
                     if (TSCAL == ONE) GO TO 135;
                  }
                     TJJ = ABS( TJJS );
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
                        X( J ) = X( J ) / TJJS;
                     } else if ( TJJ > ZERO ) {

                        // 0 < abs(A(j,j)) <= SMLNUM:

                        if ( XJ > TJJ*BIGNUM ) {

                           // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.

                           REC = ( TJJ*BIGNUM ) / XJ;
                           sscal(N, REC, X, 1 );
                           SCALE = SCALE*REC;
                           XMAX = XMAX*REC;
                        }
                        X( J ) = X( J ) / TJJS;
                     } else {

                        // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                        // scale = 0, and compute a solution to A**T*x = 0.

                        for (I = 1; I <= N; I++) { // 130
                           X( I ) = ZERO;
                        } // 130
                        X( J ) = ONE;
                        SCALE = ZERO;
                        XMAX = ZERO;
                     }
                  } // 135
               } else {

                  // Compute x(j) := x(j) / A(j,j) - sumj if the dot
                  // product has already been divided by 1/A(j,j).

                  X( J ) = X( J ) / TJJS - SUMJ;
               }
               XMAX = MAX( XMAX, ABS( X( J ) ) );
            } // 140
         }
         SCALE = SCALE / TSCAL;
      }

      // Scale the column norms by 1/TSCAL for return.

      if ( TSCAL != ONE ) {
         sscal(N, ONE / TSCAL, CNORM, 1 );
      }

      return;

      // End of SLATBS

      }
