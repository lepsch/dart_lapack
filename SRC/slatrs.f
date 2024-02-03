      SUBROUTINE SLATRS( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE, CNORM, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, NORMIN, TRANS, UPLO;
      int                INFO, LDA, N;
      REAL               SCALE
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), CNORM( * ), X( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, HALF, ONE
      const              ZERO = 0.0E+0, HALF = 0.5E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRAN, NOUNIT, UPPER;
      int                I, IMAX, J, JFIRST, JINC, JLAST;
      REAL               BIGNUM, GROW, REC, SMLNUM, SUMJ, TJJ, TJJS, TMAX, TSCAL, USCAL, XBND, XJ, XMAX
      // ..
      // .. Local Arrays ..
      REAL               WORK(1)
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ISAMAX;
      REAL               SASUM, SDOT, SLAMCH, SLANGE
      // EXTERNAL LSAME, ISAMAX, SASUM, SDOT, SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SSCAL, STRSV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOTRAN = LSAME( TRANS, 'N' )
      NOUNIT = LSAME( DIAG, 'N' )

      // Test the input parameters.

      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. LSAME( TRANS, 'C' ) ) {
         INFO = -2
      } else if ( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) {
         INFO = -3
      } else if ( .NOT.LSAME( NORMIN, 'Y' ) .AND. .NOT. LSAME( NORMIN, 'N' ) ) {
         INFO = -4
      } else if ( N.LT.0 ) {
         INFO = -5
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -7
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'SLATRS', -INFO )
         RETURN
      }

      // Quick return if possible

      SCALE = ONE
      IF( N.EQ.0 ) RETURN

      // Determine machine dependent parameters to control overflow.

      SMLNUM = SLAMCH( 'Safe minimum' ) / SLAMCH( 'Precision' )
      BIGNUM = ONE / SMLNUM

      if ( LSAME( NORMIN, 'N' ) ) {

         // Compute the 1-norm of each column, not including the diagonal.

         if ( UPPER ) {

            // A is upper triangular.

            DO 10 J = 1, N
               CNORM( J ) = SASUM( J-1, A( 1, J ), 1 )
   10       CONTINUE
         } else {

            // A is lower triangular.

            DO 20 J = 1, N - 1
               CNORM( J ) = SASUM( N-J, A( J+1, J ), 1 )
   20       CONTINUE
            CNORM( N ) = ZERO
         }
      }

      // Scale the column norms by TSCAL if the maximum element in CNORM is
      // greater than BIGNUM.

      IMAX = ISAMAX( N, CNORM, 1 )
      TMAX = CNORM( IMAX )
      if ( TMAX.LE.BIGNUM ) {
         TSCAL = ONE
      } else {

         // Avoid NaN generation if entries in CNORM exceed the
         // overflow threshold

         if ( TMAX.LE.SLAMCH('Overflow') ) {
            // Case 1: All entries in CNORM are valid floating-point numbers
            TSCAL = ONE / ( SMLNUM*TMAX )
            CALL SSCAL( N, TSCAL, CNORM, 1 )
         } else {
            // Case 2: At least one column norm of A cannot be represented
            // as floating-point number. Find the offdiagonal entry A( I, J )
            // with the largest absolute value. If this entry is not +/- Infinity,
            // use this value as TSCAL.
            TMAX = ZERO
            if ( UPPER ) {

               // A is upper triangular.

               DO J = 2, N
                  TMAX = MAX( SLANGE( 'M', J-1, 1, A( 1, J ), 1, WORK ), TMAX )
               END DO
            } else {

               // A is lower triangular.

               DO J = 1, N - 1
                  TMAX = MAX( SLANGE( 'M', N-J, 1, A( J+1, J ), 1, WORK ), TMAX )
               END DO
            }

            if ( TMAX.LE.SLAMCH('Overflow') ) {
               TSCAL = ONE / ( SMLNUM*TMAX )
               DO J = 1, N
                  if ( CNORM( J ).LE.SLAMCH('Overflow') ) {
                     CNORM( J ) = CNORM( J )*TSCAL
                  } else {
                     // Recompute the 1-norm without introducing Infinity
                     // in the summation
                     CNORM( J ) = ZERO
                     if ( UPPER ) {
                        DO I = 1, J - 1
                           CNORM( J ) = CNORM( J ) + TSCAL * ABS( A( I, J ) )
                        END DO
                     } else {
                        DO I = J + 1, N
                           CNORM( J ) = CNORM( J ) + TSCAL * ABS( A( I, J ) )
                        END DO
                     }
                  }
               END DO
            } else {
               // At least one entry of A is not a valid floating-point entry.
               // Rely on TRSV to propagate Inf and NaN.
               CALL STRSV( UPLO, TRANS, DIAG, N, A, LDA, X, 1 )
               RETURN
            }
         }
      }

      // Compute a bound on the computed solution vector to see if the
      // Level 2 BLAS routine STRSV can be used.

      J = ISAMAX( N, X, 1 )
      XMAX = ABS( X( J ) )
      XBND = XMAX
      if ( NOTRAN ) {

         // Compute the growth in A * x = b.

         if ( UPPER ) {
            JFIRST = N
            JLAST = 1
            JINC = -1
         } else {
            JFIRST = 1
            JLAST = N
            JINC = 1
         }

         if ( TSCAL.NE.ONE ) {
            GROW = ZERO
            GO TO 50
         }

         if ( NOUNIT ) {

            // A is non-unit triangular.

            // Compute GROW = 1/G(j) and XBND = 1/M(j).
            // Initially, G(0) = max{x(i), i=1,...,n}.

            GROW = ONE / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 30 J = JFIRST, JLAST, JINC

               // Exit the loop if the growth factor is too small.

               IF( GROW.LE.SMLNUM ) GO TO 50

               // M(j) = G(j-1) / abs(A(j,j))

               TJJ = ABS( A( J, J ) )
               XBND = MIN( XBND, MIN( ONE, TJJ )*GROW )
               if ( TJJ+CNORM( J ).GE.SMLNUM ) {

                  // G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )

                  GROW = GROW*( TJJ / ( TJJ+CNORM( J ) ) )
               } else {

                  // G(j) could overflow, set GROW to 0.

                  GROW = ZERO
               }
   30       CONTINUE
            GROW = XBND
         } else {

            // A is unit triangular.

            // Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.

            GROW = MIN( ONE, ONE / MAX( XBND, SMLNUM ) )
            DO 40 J = JFIRST, JLAST, JINC

               // Exit the loop if the growth factor is too small.

               IF( GROW.LE.SMLNUM ) GO TO 50

               // G(j) = G(j-1)*( 1 + CNORM(j) )

               GROW = GROW*( ONE / ( ONE+CNORM( J ) ) )
   40       CONTINUE
         }
   50    CONTINUE

      } else {

         // Compute the growth in A**T * x = b.

         if ( UPPER ) {
            JFIRST = 1
            JLAST = N
            JINC = 1
         } else {
            JFIRST = N
            JLAST = 1
            JINC = -1
         }

         if ( TSCAL.NE.ONE ) {
            GROW = ZERO
            GO TO 80
         }

         if ( NOUNIT ) {

            // A is non-unit triangular.

            // Compute GROW = 1/G(j) and XBND = 1/M(j).
            // Initially, M(0) = max{x(i), i=1,...,n}.

            GROW = ONE / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 60 J = JFIRST, JLAST, JINC

               // Exit the loop if the growth factor is too small.

               IF( GROW.LE.SMLNUM ) GO TO 80

               // G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )

               XJ = ONE + CNORM( J )
               GROW = MIN( GROW, XBND / XJ )

               // M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))

               TJJ = ABS( A( J, J ) )
               IF( XJ.GT.TJJ ) XBND = XBND*( TJJ / XJ )
   60       CONTINUE
            GROW = MIN( GROW, XBND )
         } else {

            // A is unit triangular.

            // Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.

            GROW = MIN( ONE, ONE / MAX( XBND, SMLNUM ) )
            DO 70 J = JFIRST, JLAST, JINC

               // Exit the loop if the growth factor is too small.

               IF( GROW.LE.SMLNUM ) GO TO 80

               // G(j) = ( 1 + CNORM(j) )*G(j-1)

               XJ = ONE + CNORM( J )
               GROW = GROW / XJ
   70       CONTINUE
         }
   80    CONTINUE
      }

      if ( ( GROW*TSCAL ).GT.SMLNUM ) {

         // Use the Level 2 BLAS solve if the reciprocal of the bound on
         // elements of X is not too small.

         CALL STRSV( UPLO, TRANS, DIAG, N, A, LDA, X, 1 )
      } else {

         // Use a Level 1 BLAS solve, scaling intermediate results.

         if ( XMAX.GT.BIGNUM ) {

            // Scale X so that its components are less than or equal to
            // BIGNUM in absolute value.

            SCALE = BIGNUM / XMAX
            CALL SSCAL( N, SCALE, X, 1 )
            XMAX = BIGNUM
         }

         if ( NOTRAN ) {

            // Solve A * x = b

            DO 100 J = JFIRST, JLAST, JINC

               // Compute x(j) = b(j) / A(j,j), scaling x if necessary.

               XJ = ABS( X( J ) )
               if ( NOUNIT ) {
                  TJJS = A( J, J )*TSCAL
               } else {
                  TJJS = TSCAL
                  IF( TSCAL.EQ.ONE ) GO TO 95
               }
                  TJJ = ABS( TJJS )
                  if ( TJJ.GT.SMLNUM ) {

                     // abs(A(j,j)) > SMLNUM:

                     if ( TJJ.LT.ONE ) {
                        if ( XJ.GT.TJJ*BIGNUM ) {

                           // Scale x by 1/b(j).

                           REC = ONE / XJ
                           CALL SSCAL( N, REC, X, 1 )
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        }
                     }
                     X( J ) = X( J ) / TJJS
                     XJ = ABS( X( J ) )
                  } else if ( TJJ.GT.ZERO ) {

                     // 0 < abs(A(j,j)) <= SMLNUM:

                     if ( XJ.GT.TJJ*BIGNUM ) {

                        // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
                       t // o avoid overflow when dividing by A(j,j).

                        REC = ( TJJ*BIGNUM ) / XJ
                        if ( CNORM( J ).GT.ONE ) {

                           // Scale by 1/CNORM(j) to avoid overflow when
                           // multiplying x(j) times column j.

                           REC = REC / CNORM( J )
                        }
                        CALL SSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     }
                     X( J ) = X( J ) / TJJS
                     XJ = ABS( X( J ) )
                  } else {

                     // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                     // scale = 0, and compute a solution to A*x = 0.

                     DO 90 I = 1, N
                        X( I ) = ZERO
   90                CONTINUE
                     X( J ) = ONE
                     XJ = ONE
                     SCALE = ZERO
                     XMAX = ZERO
                  }
   95          CONTINUE

               // Scale x if necessary to avoid overflow when adding a
               // multiple of column j of A.

               if ( XJ.GT.ONE ) {
                  REC = ONE / XJ
                  if ( CNORM( J ).GT.( BIGNUM-XMAX )*REC ) {

                     // Scale x by 1/(2*abs(x(j))).

                     REC = REC*HALF
                     CALL SSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                  }
               } else if ( XJ*CNORM( J ).GT.( BIGNUM-XMAX ) ) {

                  // Scale x by 1/2.

                  CALL SSCAL( N, HALF, X, 1 )
                  SCALE = SCALE*HALF
               }

               if ( UPPER ) {
                  if ( J.GT.1 ) {

                     // Compute the update
                        // x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)

                     CALL SAXPY( J-1, -X( J )*TSCAL, A( 1, J ), 1, X, 1 )
                     I = ISAMAX( J-1, X, 1 )
                     XMAX = ABS( X( I ) )
                  }
               } else {
                  if ( J.LT.N ) {

                     // Compute the update
                        // x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)

                     CALL SAXPY( N-J, -X( J )*TSCAL, A( J+1, J ), 1, X( J+1 ), 1 )
                     I = J + ISAMAX( N-J, X( J+1 ), 1 )
                     XMAX = ABS( X( I ) )
                  }
               }
  100       CONTINUE

         } else {

            // Solve A**T * x = b

            DO 140 J = JFIRST, JLAST, JINC

               // Compute x(j) = b(j) - sum A(k,j)*x(k).
                                     // k<>j

               XJ = ABS( X( J ) )
               USCAL = TSCAL
               REC = ONE / MAX( XMAX, ONE )
               if ( CNORM( J ).GT.( BIGNUM-XJ )*REC ) {

                  // If x(j) could overflow, scale x by 1/(2*XMAX).

                  REC = REC*HALF
                  if ( NOUNIT ) {
                     TJJS = A( J, J )*TSCAL
                  } else {
                     TJJS = TSCAL
                  }
                     TJJ = ABS( TJJS )
                     if ( TJJ.GT.ONE ) {

                        // Divide by A(j,j) when scaling x if A(j,j) > 1.

                        REC = MIN( ONE, REC*TJJ )
                        USCAL = USCAL / TJJS
                     }
                  if ( REC.LT.ONE ) {
                     CALL SSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  }
               }

               SUMJ = ZERO
               if ( USCAL.EQ.ONE ) {

                  // If the scaling needed for A in the dot product is 1,
                  // call SDOT to perform the dot product.

                  if ( UPPER ) {
                     SUMJ = SDOT( J-1, A( 1, J ), 1, X, 1 )
                  } else if ( J.LT.N ) {
                     SUMJ = SDOT( N-J, A( J+1, J ), 1, X( J+1 ), 1 )
                  }
               } else {

                  // Otherwise, use in-line code for the dot product.

                  if ( UPPER ) {
                     DO 110 I = 1, J - 1
                        SUMJ = SUMJ + ( A( I, J )*USCAL )*X( I )
  110                CONTINUE
                  } else if ( J.LT.N ) {
                     DO 120 I = J + 1, N
                        SUMJ = SUMJ + ( A( I, J )*USCAL )*X( I )
  120                CONTINUE
                  }
               }

               if ( USCAL.EQ.TSCAL ) {

                  // Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j)
                  // was not used to scale the dotproduct.

                  X( J ) = X( J ) - SUMJ
                  XJ = ABS( X( J ) )
                  if ( NOUNIT ) {
                     TJJS = A( J, J )*TSCAL
                  } else {
                     TJJS = TSCAL
                     IF( TSCAL.EQ.ONE ) GO TO 135
                  }

                     // Compute x(j) = x(j) / A(j,j), scaling if necessary.

                     TJJ = ABS( TJJS )
                     if ( TJJ.GT.SMLNUM ) {

                        // abs(A(j,j)) > SMLNUM:

                        if ( TJJ.LT.ONE ) {
                           if ( XJ.GT.TJJ*BIGNUM ) {

                              // Scale X by 1/abs(x(j)).

                              REC = ONE / XJ
                              CALL SSCAL( N, REC, X, 1 )
                              SCALE = SCALE*REC
                              XMAX = XMAX*REC
                           }
                        }
                        X( J ) = X( J ) / TJJS
                     } else if ( TJJ.GT.ZERO ) {

                        // 0 < abs(A(j,j)) <= SMLNUM:

                        if ( XJ.GT.TJJ*BIGNUM ) {

                           // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.

                           REC = ( TJJ*BIGNUM ) / XJ
                           CALL SSCAL( N, REC, X, 1 )
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        }
                        X( J ) = X( J ) / TJJS
                     } else {

                        // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                        // scale = 0, and compute a solution to A**T*x = 0.

                        DO 130 I = 1, N
                           X( I ) = ZERO
  130                   CONTINUE
                        X( J ) = ONE
                        SCALE = ZERO
                        XMAX = ZERO
                     }
  135             CONTINUE
               } else {

                  // Compute x(j) := x(j) / A(j,j)  - sumj if the dot
                  // product has already been divided by 1/A(j,j).

                  X( J ) = X( J ) / TJJS - SUMJ
               }
               XMAX = MAX( XMAX, ABS( X( J ) ) )
  140       CONTINUE
         }
         SCALE = SCALE / TSCAL
      }

      // Scale the column norms by 1/TSCAL for return.

      if ( TSCAL.NE.ONE ) {
         CALL SSCAL( N, ONE / TSCAL, CNORM, 1 )
      }

      RETURN

      // End of SLATRS

      }
