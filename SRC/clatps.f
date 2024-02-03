      SUBROUTINE CLATPS( UPLO, TRANS, DIAG, NORMIN, N, AP, X, SCALE, CNORM, INFO );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, NORMIN, TRANS, UPLO;
      int                INFO, N;
      REAL               SCALE;
      // ..
      // .. Array Arguments ..
      REAL               CNORM( * );
      COMPLEX            AP( * ), X( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, HALF, ONE, TWO;
      const              ZERO = 0.0, HALF = 0.5, ONE = 1.0, TWO = 2.0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRAN, NOUNIT, UPPER;
      int                I, IMAX, IP, J, JFIRST, JINC, JLAST, JLEN;
      REAL               BIGNUM, GROW, REC, SMLNUM, TJJ, TMAX, TSCAL, XBND, XJ, XMAX;
      COMPLEX            CSUMJ, TJJS, USCAL, ZDUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ICAMAX, ISAMAX;
      REAL               SCASUM, SLAMCH;
      COMPLEX            CDOTC, CDOTU, CLADIV;
      // EXTERNAL LSAME, ICAMAX, ISAMAX, SCASUM, SLAMCH, CDOTC, CDOTU, CLADIV
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CSSCAL, CTPSV, SSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CMPLX, CONJG, MAX, MIN, REAL
      // ..
      // .. Statement Functions ..
      REAL               CABS1, CABS2;
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) );
      CABS2( ZDUM ) = ABS( REAL( ZDUM ) / 2. ) + ABS( AIMAG( ZDUM ) / 2. );
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
      }
      if ( INFO != 0 ) {
         xerbla('CLATPS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Determine machine dependent parameters to control overflow.

      SMLNUM = SLAMCH( 'Safe minimum' );
      BIGNUM = ONE / SMLNUM;
      SMLNUM = SMLNUM / SLAMCH( 'Precision' );
      BIGNUM = ONE / SMLNUM;
      SCALE = ONE;

      if ( LSAME( NORMIN, 'N' ) ) {

         // Compute the 1-norm of each column, not including the diagonal.

         if ( UPPER ) {

            // A is upper triangular.

            IP = 1;
            for (J = 1; J <= N; J++) { // 10
               CNORM( J ) = SCASUM( J-1, AP( IP ), 1 );
               IP = IP + J;
            } // 10
         } else {

            // A is lower triangular.

            IP = 1;
            for (J = 1; J <= N - 1; J++) { // 20
               CNORM( J ) = SCASUM( N-J, AP( IP+1 ), 1 );
               IP = IP + N - J + 1;
            } // 20
            CNORM( N ) = ZERO;
         }
      }

      // Scale the column norms by TSCAL if the maximum element in CNORM is
      // greater than BIGNUM/2.

      IMAX = ISAMAX( N, CNORM, 1 );
      TMAX = CNORM( IMAX );
      if ( TMAX <= BIGNUM*HALF ) {
         TSCAL = ONE;
      } else {
         TSCAL = HALF / ( SMLNUM*TMAX );
         sscal(N, TSCAL, CNORM, 1 );
      }

      // Compute a bound on the computed solution vector to see if the
      // Level 2 BLAS routine CTPSV can be used.

      XMAX = ZERO;
      for (J = 1; J <= N; J++) { // 30
         XMAX = MAX( XMAX, CABS2( X( J ) ) );
      } // 30
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
            GO TO 60;
         }

         if ( NOUNIT ) {

            // A is non-unit triangular.

            // Compute GROW = 1/G(j) and XBND = 1/M(j).
            // Initially, G(0) = max{x(i), i=1,...,n}.

            GROW = HALF / MAX( XBND, SMLNUM );
            XBND = GROW;
            IP = JFIRST*( JFIRST+1 ) / 2;
            JLEN = N;
            DO 40 J = JFIRST, JLAST, JINC;

               // Exit the loop if the growth factor is too small.

               if (GROW <= SMLNUM) GO TO 60;

               TJJS = AP( IP );
               TJJ = CABS1( TJJS );

               if ( TJJ >= SMLNUM ) {

                  // M(j) = G(j-1) / abs(A(j,j))

                  XBND = MIN( XBND, MIN( ONE, TJJ )*GROW );
               } else {

                  // M(j) could overflow, set XBND to 0.

                  XBND = ZERO;
               }

               if ( TJJ+CNORM( J ) >= SMLNUM ) {

                  // G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )

                  GROW = GROW*( TJJ / ( TJJ+CNORM( J ) ) );
               } else {

                  // G(j) could overflow, set GROW to 0.

                  GROW = ZERO;
               }
               IP = IP + JINC*JLEN;
               JLEN = JLEN - 1;
            } // 40
            GROW = XBND;
         } else {

            // A is unit triangular.

            // Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.

            GROW = MIN( ONE, HALF / MAX( XBND, SMLNUM ) );
            DO 50 J = JFIRST, JLAST, JINC;

               // Exit the loop if the growth factor is too small.

               if (GROW <= SMLNUM) GO TO 60;

               // G(j) = G(j-1)*( 1 + CNORM(j) )

               GROW = GROW*( ONE / ( ONE+CNORM( J ) ) );
            } // 50
         }
         } // 60

      } else {

         // Compute the growth in A**T * x = b  or  A**H * x = b.

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
            GO TO 90;
         }

         if ( NOUNIT ) {

            // A is non-unit triangular.

            // Compute GROW = 1/G(j) and XBND = 1/M(j).
            // Initially, M(0) = max{x(i), i=1,...,n}.

            GROW = HALF / MAX( XBND, SMLNUM );
            XBND = GROW;
            IP = JFIRST*( JFIRST+1 ) / 2;
            JLEN = 1;
            DO 70 J = JFIRST, JLAST, JINC;

               // Exit the loop if the growth factor is too small.

               if (GROW <= SMLNUM) GO TO 90;

               // G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )

               XJ = ONE + CNORM( J );
               GROW = MIN( GROW, XBND / XJ );

               TJJS = AP( IP );
               TJJ = CABS1( TJJS );

               if ( TJJ >= SMLNUM ) {

                  // M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))

                  if (XJ > TJJ) XBND = XBND*( TJJ / XJ );
               } else {

                  // M(j) could overflow, set XBND to 0.

                  XBND = ZERO;
               }
               JLEN = JLEN + 1;
               IP = IP + JINC*JLEN;
            } // 70
            GROW = MIN( GROW, XBND );
         } else {

            // A is unit triangular.

            // Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.

            GROW = MIN( ONE, HALF / MAX( XBND, SMLNUM ) );
            DO 80 J = JFIRST, JLAST, JINC;

               // Exit the loop if the growth factor is too small.

               if (GROW <= SMLNUM) GO TO 90;

               // G(j) = ( 1 + CNORM(j) )*G(j-1)

               XJ = ONE + CNORM( J );
               GROW = GROW / XJ;
            } // 80
         }
         } // 90
      }

      if ( ( GROW*TSCAL ) > SMLNUM ) {

         // Use the Level 2 BLAS solve if the reciprocal of the bound on
         // elements of X is not too small.

         ctpsv(UPLO, TRANS, DIAG, N, AP, X, 1 );
      } else {

         // Use a Level 1 BLAS solve, scaling intermediate results.

         if ( XMAX > BIGNUM*HALF ) {

            // Scale X so that its components are less than or equal to
            // BIGNUM in absolute value.

            SCALE = ( BIGNUM*HALF ) / XMAX;
            csscal(N, SCALE, X, 1 );
            XMAX = BIGNUM;
         } else {
            XMAX = XMAX*TWO;
         }

         if ( NOTRAN ) {

            // Solve A * x = b

            IP = JFIRST*( JFIRST+1 ) / 2;
            DO 110 J = JFIRST, JLAST, JINC;

               // Compute x(j) = b(j) / A(j,j), scaling x if necessary.

               XJ = CABS1( X( J ) );
               if ( NOUNIT ) {
                  TJJS = AP( IP )*TSCAL;
               } else {
                  TJJS = TSCAL;
                  if (TSCAL == ONE) GO TO 105;
               }
                  TJJ = CABS1( TJJS );
                  if ( TJJ > SMLNUM ) {

                     // abs(A(j,j)) > SMLNUM:

                     if ( TJJ < ONE ) {
                        if ( XJ > TJJ*BIGNUM ) {

                           // Scale x by 1/b(j).

                           REC = ONE / XJ;
                           csscal(N, REC, X, 1 );
                           SCALE = SCALE*REC;
                           XMAX = XMAX*REC;
                        }
                     }
                     X( J ) = CLADIV( X( J ), TJJS );
                     XJ = CABS1( X( J ) );
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
                        csscal(N, REC, X, 1 );
                        SCALE = SCALE*REC;
                        XMAX = XMAX*REC;
                     }
                     X( J ) = CLADIV( X( J ), TJJS );
                     XJ = CABS1( X( J ) );
                  } else {

                     // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                     // scale = 0, and compute a solution to A*x = 0.

                     for (I = 1; I <= N; I++) { // 100
                        X( I ) = ZERO;
                     } // 100
                     X( J ) = ONE;
                     XJ = ONE;
                     SCALE = ZERO;
                     XMAX = ZERO;
                  }
               } // 105

               // Scale x if necessary to avoid overflow when adding a
               // multiple of column j of A.

               if ( XJ > ONE ) {
                  REC = ONE / XJ;
                  if ( CNORM( J ) > ( BIGNUM-XMAX )*REC ) {

                     // Scale x by 1/(2*abs(x(j))).

                     REC = REC*HALF;
                     csscal(N, REC, X, 1 );
                     SCALE = SCALE*REC;
                  }
               } else if ( XJ*CNORM( J ) > ( BIGNUM-XMAX ) ) {

                  // Scale x by 1/2.

                  csscal(N, HALF, X, 1 );
                  SCALE = SCALE*HALF;
               }

               if ( UPPER ) {
                  if ( J > 1 ) {

                     // Compute the update
                        // x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)

                     caxpy(J-1, -X( J )*TSCAL, AP( IP-J+1 ), 1, X, 1 );
                     I = ICAMAX( J-1, X, 1 );
                     XMAX = CABS1( X( I ) );
                  }
                  IP = IP - J;
               } else {
                  if ( J < N ) {

                     // Compute the update
                        // x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)

                     caxpy(N-J, -X( J )*TSCAL, AP( IP+1 ), 1, X( J+1 ), 1 );
                     I = J + ICAMAX( N-J, X( J+1 ), 1 );
                     XMAX = CABS1( X( I ) );
                  }
                  IP = IP + N - J + 1;
               }
            } // 110

         } else if ( LSAME( TRANS, 'T' ) ) {

            // Solve A**T * x = b

            IP = JFIRST*( JFIRST+1 ) / 2;
            JLEN = 1;
            DO 150 J = JFIRST, JLAST, JINC;

               // Compute x(j) = b(j) - sum A(k,j)*x(k).
                                     // k<>j

               XJ = CABS1( X( J ) );
               USCAL = TSCAL;
               REC = ONE / MAX( XMAX, ONE );
               if ( CNORM( J ) > ( BIGNUM-XJ )*REC ) {

                  // If x(j) could overflow, scale x by 1/(2*XMAX).

                  REC = REC*HALF;
                  if ( NOUNIT ) {
                     TJJS = AP( IP )*TSCAL;
                  } else {
                     TJJS = TSCAL;
                  }
                     TJJ = CABS1( TJJS );
                     if ( TJJ > ONE ) {

                        // Divide by A(j,j) when scaling x if A(j,j) > 1.

                        REC = MIN( ONE, REC*TJJ );
                        USCAL = CLADIV( USCAL, TJJS );
                     }
                  if ( REC < ONE ) {
                     csscal(N, REC, X, 1 );
                     SCALE = SCALE*REC;
                     XMAX = XMAX*REC;
                  }
               }

               CSUMJ = ZERO;
               if ( USCAL == CMPLX( ONE ) ) {

                  // If the scaling needed for A in the dot product is 1,
                  // call CDOTU to perform the dot product.

                  if ( UPPER ) {
                     CSUMJ = CDOTU( J-1, AP( IP-J+1 ), 1, X, 1 );
                  } else if ( J < N ) {
                     CSUMJ = CDOTU( N-J, AP( IP+1 ), 1, X( J+1 ), 1 );
                  }
               } else {

                  // Otherwise, use in-line code for the dot product.

                  if ( UPPER ) {
                     for (I = 1; I <= J - 1; I++) { // 120
                        CSUMJ = CSUMJ + ( AP( IP-J+I )*USCAL )*X( I );
                     } // 120
                  } else if ( J < N ) {
                     for (I = 1; I <= N - J; I++) { // 130
                        CSUMJ = CSUMJ + ( AP( IP+I )*USCAL )*X( J+I );
                     } // 130
                  }
               }

               if ( USCAL == CMPLX( TSCAL ) ) {

                  // Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
                  // was not used to scale the dotproduct.

                  X( J ) = X( J ) - CSUMJ;
                  XJ = CABS1( X( J ) );
                  if ( NOUNIT ) {

                     // Compute x(j) = x(j) / A(j,j), scaling if necessary.

                     TJJS = AP( IP )*TSCAL;
                  } else {
                     TJJS = TSCAL;
                     if (TSCAL == ONE) GO TO 145;
                  }
                     TJJ = CABS1( TJJS );
                     if ( TJJ > SMLNUM ) {

                        // abs(A(j,j)) > SMLNUM:

                        if ( TJJ < ONE ) {
                           if ( XJ > TJJ*BIGNUM ) {

                              // Scale X by 1/abs(x(j)).

                              REC = ONE / XJ;
                              csscal(N, REC, X, 1 );
                              SCALE = SCALE*REC;
                              XMAX = XMAX*REC;
                           }
                        }
                        X( J ) = CLADIV( X( J ), TJJS );
                     } else if ( TJJ > ZERO ) {

                        // 0 < abs(A(j,j)) <= SMLNUM:

                        if ( XJ > TJJ*BIGNUM ) {

                           // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.

                           REC = ( TJJ*BIGNUM ) / XJ;
                           csscal(N, REC, X, 1 );
                           SCALE = SCALE*REC;
                           XMAX = XMAX*REC;
                        }
                        X( J ) = CLADIV( X( J ), TJJS );
                     } else {

                        // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                        // scale = 0 and compute a solution to A**T *x = 0.

                        for (I = 1; I <= N; I++) { // 140
                           X( I ) = ZERO;
                        } // 140
                        X( J ) = ONE;
                        SCALE = ZERO;
                        XMAX = ZERO;
                     }
                  } // 145
               } else {

                  // Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
                  // product has already been divided by 1/A(j,j).

                  X( J ) = CLADIV( X( J ), TJJS ) - CSUMJ;
               }
               XMAX = MAX( XMAX, CABS1( X( J ) ) );
               JLEN = JLEN + 1;
               IP = IP + JINC*JLEN;
            } // 150

         } else {

            // Solve A**H * x = b

            IP = JFIRST*( JFIRST+1 ) / 2;
            JLEN = 1;
            DO 190 J = JFIRST, JLAST, JINC;

               // Compute x(j) = b(j) - sum A(k,j)*x(k).
                                     // k<>j

               XJ = CABS1( X( J ) );
               USCAL = TSCAL;
               REC = ONE / MAX( XMAX, ONE );
               if ( CNORM( J ) > ( BIGNUM-XJ )*REC ) {

                  // If x(j) could overflow, scale x by 1/(2*XMAX).

                  REC = REC*HALF;
                  if ( NOUNIT ) {
                     TJJS = CONJG( AP( IP ) )*TSCAL;
                  } else {
                     TJJS = TSCAL;
                  }
                     TJJ = CABS1( TJJS );
                     if ( TJJ > ONE ) {

                        // Divide by A(j,j) when scaling x if A(j,j) > 1.

                        REC = MIN( ONE, REC*TJJ );
                        USCAL = CLADIV( USCAL, TJJS );
                     }
                  if ( REC < ONE ) {
                     csscal(N, REC, X, 1 );
                     SCALE = SCALE*REC;
                     XMAX = XMAX*REC;
                  }
               }

               CSUMJ = ZERO;
               if ( USCAL == CMPLX( ONE ) ) {

                  // If the scaling needed for A in the dot product is 1,
                  // call CDOTC to perform the dot product.

                  if ( UPPER ) {
                     CSUMJ = CDOTC( J-1, AP( IP-J+1 ), 1, X, 1 );
                  } else if ( J < N ) {
                     CSUMJ = CDOTC( N-J, AP( IP+1 ), 1, X( J+1 ), 1 );
                  }
               } else {

                  // Otherwise, use in-line code for the dot product.

                  if ( UPPER ) {
                     for (I = 1; I <= J - 1; I++) { // 160
                        CSUMJ = CSUMJ + ( CONJG( AP( IP-J+I ) )*USCAL )* X( I );
                     } // 160
                  } else if ( J < N ) {
                     for (I = 1; I <= N - J; I++) { // 170
                        CSUMJ = CSUMJ + ( CONJG( AP( IP+I ) )*USCAL )* X( J+I );
                     } // 170
                  }
               }

               if ( USCAL == CMPLX( TSCAL ) ) {

                  // Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
                  // was not used to scale the dotproduct.

                  X( J ) = X( J ) - CSUMJ;
                  XJ = CABS1( X( J ) );
                  if ( NOUNIT ) {

                     // Compute x(j) = x(j) / A(j,j), scaling if necessary.

                     TJJS = CONJG( AP( IP ) )*TSCAL;
                  } else {
                     TJJS = TSCAL;
                     if (TSCAL == ONE) GO TO 185;
                  }
                     TJJ = CABS1( TJJS );
                     if ( TJJ > SMLNUM ) {

                        // abs(A(j,j)) > SMLNUM:

                        if ( TJJ < ONE ) {
                           if ( XJ > TJJ*BIGNUM ) {

                              // Scale X by 1/abs(x(j)).

                              REC = ONE / XJ;
                              csscal(N, REC, X, 1 );
                              SCALE = SCALE*REC;
                              XMAX = XMAX*REC;
                           }
                        }
                        X( J ) = CLADIV( X( J ), TJJS );
                     } else if ( TJJ > ZERO ) {

                        // 0 < abs(A(j,j)) <= SMLNUM:

                        if ( XJ > TJJ*BIGNUM ) {

                           // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.

                           REC = ( TJJ*BIGNUM ) / XJ;
                           csscal(N, REC, X, 1 );
                           SCALE = SCALE*REC;
                           XMAX = XMAX*REC;
                        }
                        X( J ) = CLADIV( X( J ), TJJS );
                     } else {

                        // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                        // scale = 0 and compute a solution to A**H *x = 0.

                        for (I = 1; I <= N; I++) { // 180
                           X( I ) = ZERO;
                        } // 180
                        X( J ) = ONE;
                        SCALE = ZERO;
                        XMAX = ZERO;
                     }
                  } // 185
               } else {

                  // Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
                  // product has already been divided by 1/A(j,j).

                  X( J ) = CLADIV( X( J ), TJJS ) - CSUMJ;
               }
               XMAX = MAX( XMAX, CABS1( X( J ) ) );
               JLEN = JLEN + 1;
               IP = IP + JINC*JLEN;
            } // 190
         }
         SCALE = SCALE / TSCAL;
      }

      // Scale the column norms by 1/TSCAL for return.

      if ( TSCAL != ONE ) {
         sscal(N, ONE / TSCAL, CNORM, 1 );
      }

      return;

      // End of CLATPS

      }
