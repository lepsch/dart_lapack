      SUBROUTINE CLATBS( UPLO, TRANS, DIAG, NORMIN, N, KD, AB, LDAB, X, SCALE, CNORM, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, NORMIN, TRANS, UPLO;
      int                INFO, KD, LDAB, N;
      REAL               SCALE
      // ..
      // .. Array Arguments ..
      REAL               CNORM( * )
      COMPLEX            AB( LDAB, * ), X( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, HALF, ONE, TWO
      const              ZERO = 0.0E+0, HALF = 0.5E+0, ONE = 1.0E+0, TWO = 2.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRAN, NOUNIT, UPPER;
      int                I, IMAX, J, JFIRST, JINC, JLAST, JLEN, MAIND;
      REAL               BIGNUM, GROW, REC, SMLNUM, TJJ, TMAX, TSCAL, XBND, XJ, XMAX
      COMPLEX            CSUMJ, TJJS, USCAL, ZDUM
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ICAMAX, ISAMAX;
      REAL               SCASUM, SLAMCH
      COMPLEX            CDOTC, CDOTU, CLADIV
      // EXTERNAL LSAME, ICAMAX, ISAMAX, SCASUM, SLAMCH, CDOTC, CDOTU, CLADIV
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CSSCAL, CTBSV, SSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CMPLX, CONJG, MAX, MIN, REAL
      // ..
      // .. Statement Functions ..
      REAL               CABS1, CABS2
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
      CABS2( ZDUM ) = ABS( REAL( ZDUM ) / 2. ) + ABS( AIMAG( ZDUM ) / 2. )
      // ..
      // .. Executable Statements ..

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOTRAN = LSAME( TRANS, 'N' )
      NOUNIT = LSAME( DIAG, 'N' )

      // Test the input parameters.

      if ( .NOT.UPPER && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( .NOT.NOTRAN && .NOT.LSAME( TRANS, 'T' ) && .NOT. LSAME( TRANS, 'C' ) ) {
         INFO = -2
      } else if ( .NOT.NOUNIT && .NOT.LSAME( DIAG, 'U' ) ) {
         INFO = -3
      } else if ( .NOT.LSAME( NORMIN, 'Y' ) && .NOT. LSAME( NORMIN, 'N' ) ) {
         INFO = -4
      } else if ( N.LT.0 ) {
         INFO = -5
      } else if ( KD.LT.0 ) {
         INFO = -6
      } else if ( LDAB.LT.KD+1 ) {
         INFO = -8
      }
      if ( INFO != 0 ) {
         xerbla('CLATBS', -INFO );
         RETURN
      }

      // Quick return if possible

      SCALE = ONE
      if (N == 0) RETURN;

      // Determine machine dependent parameters to control overflow.

      SMLNUM = SLAMCH( 'Safe minimum' ) / SLAMCH( 'Precision' )
      BIGNUM = ONE / SMLNUM

      if ( LSAME( NORMIN, 'N' ) ) {

         // Compute the 1-norm of each column, not including the diagonal.

         if ( UPPER ) {

            // A is upper triangular.

            for (J = 1; J <= N; J++) { // 10
               JLEN = MIN( KD, J-1 )
               CNORM( J ) = SCASUM( JLEN, AB( KD+1-JLEN, J ), 1 )
            } // 10
         } else {

            // A is lower triangular.

            for (J = 1; J <= N; J++) { // 20
               JLEN = MIN( KD, N-J )
               if ( JLEN.GT.0 ) {
                  CNORM( J ) = SCASUM( JLEN, AB( 2, J ), 1 )
               } else {
                  CNORM( J ) = ZERO
               }
            } // 20
         }
      }

      // Scale the column norms by TSCAL if the maximum element in CNORM is
      // greater than BIGNUM/2.

      IMAX = ISAMAX( N, CNORM, 1 )
      TMAX = CNORM( IMAX )
      if ( TMAX.LE.BIGNUM*HALF ) {
         TSCAL = ONE
      } else {
         TSCAL = HALF / ( SMLNUM*TMAX )
         sscal(N, TSCAL, CNORM, 1 );
      }

      // Compute a bound on the computed solution vector to see if the
      // Level 2 BLAS routine CTBSV can be used.

      XMAX = ZERO
      for (J = 1; J <= N; J++) { // 30
         XMAX = MAX( XMAX, CABS2( X( J ) ) )
      } // 30
      XBND = XMAX
      if ( NOTRAN ) {

         // Compute the growth in A * x = b.

         if ( UPPER ) {
            JFIRST = N
            JLAST = 1
            JINC = -1
            MAIND = KD + 1
         } else {
            JFIRST = 1
            JLAST = N
            JINC = 1
            MAIND = 1
         }

         if ( TSCAL != ONE ) {
            GROW = ZERO
            GO TO 60
         }

         if ( NOUNIT ) {

            // A is non-unit triangular.

            // Compute GROW = 1/G(j) and XBND = 1/M(j).
            // Initially, G(0) = max{x(i), i=1,...,n}.

            GROW = HALF / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 40 J = JFIRST, JLAST, JINC

               // Exit the loop if the growth factor is too small.

               if (GROW.LE.SMLNUM) GO TO 60;

               TJJS = AB( MAIND, J )
               TJJ = CABS1( TJJS )

               if ( TJJ.GE.SMLNUM ) {

                  // M(j) = G(j-1) / abs(A(j,j))

                  XBND = MIN( XBND, MIN( ONE, TJJ )*GROW )
               } else {

                  // M(j) could overflow, set XBND to 0.

                  XBND = ZERO
               }

               if ( TJJ+CNORM( J ).GE.SMLNUM ) {

                  // G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )

                  GROW = GROW*( TJJ / ( TJJ+CNORM( J ) ) )
               } else {

                  // G(j) could overflow, set GROW to 0.

                  GROW = ZERO
               }
            } // 40
            GROW = XBND
         } else {

            // A is unit triangular.

            // Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.

            GROW = MIN( ONE, HALF / MAX( XBND, SMLNUM ) )
            DO 50 J = JFIRST, JLAST, JINC

               // Exit the loop if the growth factor is too small.

               if (GROW.LE.SMLNUM) GO TO 60;

               // G(j) = G(j-1)*( 1 + CNORM(j) )

               GROW = GROW*( ONE / ( ONE+CNORM( J ) ) )
            } // 50
         }
         } // 60

      } else {

         // Compute the growth in A**T * x = b  or  A**H * x = b.

         if ( UPPER ) {
            JFIRST = 1
            JLAST = N
            JINC = 1
            MAIND = KD + 1
         } else {
            JFIRST = N
            JLAST = 1
            JINC = -1
            MAIND = 1
         }

         if ( TSCAL != ONE ) {
            GROW = ZERO
            GO TO 90
         }

         if ( NOUNIT ) {

            // A is non-unit triangular.

            // Compute GROW = 1/G(j) and XBND = 1/M(j).
            // Initially, M(0) = max{x(i), i=1,...,n}.

            GROW = HALF / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 70 J = JFIRST, JLAST, JINC

               // Exit the loop if the growth factor is too small.

               if (GROW.LE.SMLNUM) GO TO 90;

               // G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )

               XJ = ONE + CNORM( J )
               GROW = MIN( GROW, XBND / XJ )

               TJJS = AB( MAIND, J )
               TJJ = CABS1( TJJS )

               if ( TJJ.GE.SMLNUM ) {

                  // M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))

                  if (XJ.GT.TJJ) XBND = XBND*( TJJ / XJ );
               } else {

                  // M(j) could overflow, set XBND to 0.

                  XBND = ZERO
               }
            } // 70
            GROW = MIN( GROW, XBND )
         } else {

            // A is unit triangular.

            // Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.

            GROW = MIN( ONE, HALF / MAX( XBND, SMLNUM ) )
            DO 80 J = JFIRST, JLAST, JINC

               // Exit the loop if the growth factor is too small.

               if (GROW.LE.SMLNUM) GO TO 90;

               // G(j) = ( 1 + CNORM(j) )*G(j-1)

               XJ = ONE + CNORM( J )
               GROW = GROW / XJ
            } // 80
         }
         } // 90
      }

      if ( ( GROW*TSCAL ).GT.SMLNUM ) {

         // Use the Level 2 BLAS solve if the reciprocal of the bound on
         // elements of X is not too small.

         ctbsv(UPLO, TRANS, DIAG, N, KD, AB, LDAB, X, 1 );
      } else {

         // Use a Level 1 BLAS solve, scaling intermediate results.

         if ( XMAX.GT.BIGNUM*HALF ) {

            // Scale X so that its components are less than or equal to
            // BIGNUM in absolute value.

            SCALE = ( BIGNUM*HALF ) / XMAX
            csscal(N, SCALE, X, 1 );
            XMAX = BIGNUM
         } else {
            XMAX = XMAX*TWO
         }

         if ( NOTRAN ) {

            // Solve A * x = b

            DO 110 J = JFIRST, JLAST, JINC

               // Compute x(j) = b(j) / A(j,j), scaling x if necessary.

               XJ = CABS1( X( J ) )
               if ( NOUNIT ) {
                  TJJS = AB( MAIND, J )*TSCAL
               } else {
                  TJJS = TSCAL
                  if (TSCAL == ONE) GO TO 105;
               }
                  TJJ = CABS1( TJJS )
                  if ( TJJ.GT.SMLNUM ) {

                     // abs(A(j,j)) > SMLNUM:

                     if ( TJJ.LT.ONE ) {
                        if ( XJ.GT.TJJ*BIGNUM ) {

                           // Scale x by 1/b(j).

                           REC = ONE / XJ
                           csscal(N, REC, X, 1 );
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        }
                     }
                     X( J ) = CLADIV( X( J ), TJJS )
                     XJ = CABS1( X( J ) )
                  } else if ( TJJ.GT.ZERO ) {

                     // 0 < abs(A(j,j)) <= SMLNUM:

                     if ( XJ.GT.TJJ*BIGNUM ) {

                        // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
                        // to avoid overflow when dividing by A(j,j).

                        REC = ( TJJ*BIGNUM ) / XJ
                        if ( CNORM( J ).GT.ONE ) {

                           // Scale by 1/CNORM(j) to avoid overflow when
                           // multiplying x(j) times column j.

                           REC = REC / CNORM( J )
                        }
                        csscal(N, REC, X, 1 );
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     }
                     X( J ) = CLADIV( X( J ), TJJS )
                     XJ = CABS1( X( J ) )
                  } else {

                     // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                     // scale = 0, and compute a solution to A*x = 0.

                     for (I = 1; I <= N; I++) { // 100
                        X( I ) = ZERO
                     } // 100
                     X( J ) = ONE
                     XJ = ONE
                     SCALE = ZERO
                     XMAX = ZERO
                  }
               } // 105

               // Scale x if necessary to avoid overflow when adding a
               // multiple of column j of A.

               if ( XJ.GT.ONE ) {
                  REC = ONE / XJ
                  if ( CNORM( J ).GT.( BIGNUM-XMAX )*REC ) {

                     // Scale x by 1/(2*abs(x(j))).

                     REC = REC*HALF
                     csscal(N, REC, X, 1 );
                     SCALE = SCALE*REC
                  }
               } else if ( XJ*CNORM( J ).GT.( BIGNUM-XMAX ) ) {

                  // Scale x by 1/2.

                  csscal(N, HALF, X, 1 );
                  SCALE = SCALE*HALF
               }

               if ( UPPER ) {
                  if ( J.GT.1 ) {

                     // Compute the update
                        // x(max(1,j-kd):j-1) := x(max(1,j-kd):j-1) -
                                              // x(j)* A(max(1,j-kd):j-1,j)

                     JLEN = MIN( KD, J-1 )
                     caxpy(JLEN, -X( J )*TSCAL, AB( KD+1-JLEN, J ), 1, X( J-JLEN ), 1 );
                     I = ICAMAX( J-1, X, 1 )
                     XMAX = CABS1( X( I ) )
                  }
               } else if ( J.LT.N ) {

                  // Compute the update
                     // x(j+1:min(j+kd,n)) := x(j+1:min(j+kd,n)) -
                                           // x(j) * A(j+1:min(j+kd,n),j)

                  JLEN = MIN( KD, N-J )
                  if (JLEN.GT.0) CALL CAXPY( JLEN, -X( J )*TSCAL, AB( 2, J ), 1, X( J+1 ), 1 );
                  I = J + ICAMAX( N-J, X( J+1 ), 1 )
                  XMAX = CABS1( X( I ) )
               }
            } // 110

         } else if ( LSAME( TRANS, 'T' ) ) {

            // Solve A**T * x = b

            DO 150 J = JFIRST, JLAST, JINC

               // Compute x(j) = b(j) - sum A(k,j)*x(k).
                                     // k<>j

               XJ = CABS1( X( J ) )
               USCAL = TSCAL
               REC = ONE / MAX( XMAX, ONE )
               if ( CNORM( J ).GT.( BIGNUM-XJ )*REC ) {

                  // If x(j) could overflow, scale x by 1/(2*XMAX).

                  REC = REC*HALF
                  if ( NOUNIT ) {
                     TJJS = AB( MAIND, J )*TSCAL
                  } else {
                     TJJS = TSCAL
                  }
                     TJJ = CABS1( TJJS )
                     if ( TJJ.GT.ONE ) {

                        // Divide by A(j,j) when scaling x if A(j,j) > 1.

                        REC = MIN( ONE, REC*TJJ )
                        USCAL = CLADIV( USCAL, TJJS )
                     }
                  if ( REC.LT.ONE ) {
                     csscal(N, REC, X, 1 );
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  }
               }

               CSUMJ = ZERO
               if ( USCAL == CMPLX( ONE ) ) {

                  // If the scaling needed for A in the dot product is 1,
                  // call CDOTU to perform the dot product.

                  if ( UPPER ) {
                     JLEN = MIN( KD, J-1 )
                     CSUMJ = CDOTU( JLEN, AB( KD+1-JLEN, J ), 1, X( J-JLEN ), 1 )
                  } else {
                     JLEN = MIN( KD, N-J )
                     if (JLEN.GT.1) CSUMJ = CDOTU( JLEN, AB( 2, J ), 1, X( J+1 ), 1 );
                  }
               } else {

                  // Otherwise, use in-line code for the dot product.

                  if ( UPPER ) {
                     JLEN = MIN( KD, J-1 )
                     for (I = 1; I <= JLEN; I++) { // 120
                        CSUMJ = CSUMJ + ( AB( KD+I-JLEN, J )*USCAL )* X( J-JLEN-1+I )
                     } // 120
                  } else {
                     JLEN = MIN( KD, N-J )
                     for (I = 1; I <= JLEN; I++) { // 130
                        CSUMJ = CSUMJ + ( AB( I+1, J )*USCAL )*X( J+I )
                     } // 130
                  }
               }

               if ( USCAL == CMPLX( TSCAL ) ) {

                  // Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
                  // was not used to scale the dotproduct.

                  X( J ) = X( J ) - CSUMJ
                  XJ = CABS1( X( J ) )
                  if ( NOUNIT ) {

                     // Compute x(j) = x(j) / A(j,j), scaling if necessary.

                     TJJS = AB( MAIND, J )*TSCAL
                  } else {
                     TJJS = TSCAL
                     if (TSCAL == ONE) GO TO 145;
                  }
                     TJJ = CABS1( TJJS )
                     if ( TJJ.GT.SMLNUM ) {

                        // abs(A(j,j)) > SMLNUM:

                        if ( TJJ.LT.ONE ) {
                           if ( XJ.GT.TJJ*BIGNUM ) {

                              // Scale X by 1/abs(x(j)).

                              REC = ONE / XJ
                              csscal(N, REC, X, 1 );
                              SCALE = SCALE*REC
                              XMAX = XMAX*REC
                           }
                        }
                        X( J ) = CLADIV( X( J ), TJJS )
                     } else if ( TJJ.GT.ZERO ) {

                        // 0 < abs(A(j,j)) <= SMLNUM:

                        if ( XJ.GT.TJJ*BIGNUM ) {

                           // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.

                           REC = ( TJJ*BIGNUM ) / XJ
                           csscal(N, REC, X, 1 );
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        }
                        X( J ) = CLADIV( X( J ), TJJS )
                     } else {

                        // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                        // scale = 0 and compute a solution to A**T *x = 0.

                        for (I = 1; I <= N; I++) { // 140
                           X( I ) = ZERO
                        } // 140
                        X( J ) = ONE
                        SCALE = ZERO
                        XMAX = ZERO
                     }
                  } // 145
               } else {

                  // Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
                  // product has already been divided by 1/A(j,j).

                  X( J ) = CLADIV( X( J ), TJJS ) - CSUMJ
               }
               XMAX = MAX( XMAX, CABS1( X( J ) ) )
            } // 150

         } else {

            // Solve A**H * x = b

            DO 190 J = JFIRST, JLAST, JINC

               // Compute x(j) = b(j) - sum A(k,j)*x(k).
                                     // k<>j

               XJ = CABS1( X( J ) )
               USCAL = TSCAL
               REC = ONE / MAX( XMAX, ONE )
               if ( CNORM( J ).GT.( BIGNUM-XJ )*REC ) {

                  // If x(j) could overflow, scale x by 1/(2*XMAX).

                  REC = REC*HALF
                  if ( NOUNIT ) {
                     TJJS = CONJG( AB( MAIND, J ) )*TSCAL
                  } else {
                     TJJS = TSCAL
                  }
                     TJJ = CABS1( TJJS )
                     if ( TJJ.GT.ONE ) {

                        // Divide by A(j,j) when scaling x if A(j,j) > 1.

                        REC = MIN( ONE, REC*TJJ )
                        USCAL = CLADIV( USCAL, TJJS )
                     }
                  if ( REC.LT.ONE ) {
                     csscal(N, REC, X, 1 );
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  }
               }

               CSUMJ = ZERO
               if ( USCAL == CMPLX( ONE ) ) {

                  // If the scaling needed for A in the dot product is 1,
                  // call CDOTC to perform the dot product.

                  if ( UPPER ) {
                     JLEN = MIN( KD, J-1 )
                     CSUMJ = CDOTC( JLEN, AB( KD+1-JLEN, J ), 1, X( J-JLEN ), 1 )
                  } else {
                     JLEN = MIN( KD, N-J )
                     if (JLEN.GT.1) CSUMJ = CDOTC( JLEN, AB( 2, J ), 1, X( J+1 ), 1 );
                  }
               } else {

                  // Otherwise, use in-line code for the dot product.

                  if ( UPPER ) {
                     JLEN = MIN( KD, J-1 )
                     for (I = 1; I <= JLEN; I++) { // 160
                        CSUMJ = CSUMJ + ( CONJG( AB( KD+I-JLEN, J ) )* USCAL )*X( J-JLEN-1+I )
                     } // 160
                  } else {
                     JLEN = MIN( KD, N-J )
                     for (I = 1; I <= JLEN; I++) { // 170
                        CSUMJ = CSUMJ + ( CONJG( AB( I+1, J ) )*USCAL )* X( J+I )
                     } // 170
                  }
               }

               if ( USCAL == CMPLX( TSCAL ) ) {

                  // Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
                  // was not used to scale the dotproduct.

                  X( J ) = X( J ) - CSUMJ
                  XJ = CABS1( X( J ) )
                  if ( NOUNIT ) {

                     // Compute x(j) = x(j) / A(j,j), scaling if necessary.

                     TJJS = CONJG( AB( MAIND, J ) )*TSCAL
                  } else {
                     TJJS = TSCAL
                     if (TSCAL == ONE) GO TO 185;
                  }
                     TJJ = CABS1( TJJS )
                     if ( TJJ.GT.SMLNUM ) {

                        // abs(A(j,j)) > SMLNUM:

                        if ( TJJ.LT.ONE ) {
                           if ( XJ.GT.TJJ*BIGNUM ) {

                              // Scale X by 1/abs(x(j)).

                              REC = ONE / XJ
                              csscal(N, REC, X, 1 );
                              SCALE = SCALE*REC
                              XMAX = XMAX*REC
                           }
                        }
                        X( J ) = CLADIV( X( J ), TJJS )
                     } else if ( TJJ.GT.ZERO ) {

                        // 0 < abs(A(j,j)) <= SMLNUM:

                        if ( XJ.GT.TJJ*BIGNUM ) {

                           // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.

                           REC = ( TJJ*BIGNUM ) / XJ
                           csscal(N, REC, X, 1 );
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        }
                        X( J ) = CLADIV( X( J ), TJJS )
                     } else {

                        // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                        // scale = 0 and compute a solution to A**H *x = 0.

                        for (I = 1; I <= N; I++) { // 180
                           X( I ) = ZERO
                        } // 180
                        X( J ) = ONE
                        SCALE = ZERO
                        XMAX = ZERO
                     }
                  } // 185
               } else {

                  // Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
                  // product has already been divided by 1/A(j,j).

                  X( J ) = CLADIV( X( J ), TJJS ) - CSUMJ
               }
               XMAX = MAX( XMAX, CABS1( X( J ) ) )
            } // 190
         }
         SCALE = SCALE / TSCAL
      }

      // Scale the column norms by 1/TSCAL for return.

      if ( TSCAL != ONE ) {
         sscal(N, ONE / TSCAL, CNORM, 1 );
      }

      RETURN

      // End of CLATBS

      }
