      SUBROUTINE ZLATRS( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE, CNORM, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, NORMIN, TRANS, UPLO;
      int                INFO, LDA, N;
      double             SCALE;
      // ..
      // .. Array Arguments ..
      double             CNORM( * );
      COMPLEX*16         A( LDA, * ), X( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, HALF, ONE, TWO;
      const              ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0, TWO = 2.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRAN, NOUNIT, UPPER;
      int                I, IMAX, J, JFIRST, JINC, JLAST;
      double             BIGNUM, GROW, REC, SMLNUM, TJJ, TMAX, TSCAL, XBND, XJ, XMAX;
      COMPLEX*16         CSUMJ, TJJS, USCAL, ZDUM
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IDAMAX, IZAMAX;
      double             DLAMCH, DZASUM;
      COMPLEX*16         ZDOTC, ZDOTU, ZLADIV
      // EXTERNAL LSAME, IDAMAX, IZAMAX, DLAMCH, DZASUM, ZDOTC, ZDOTU, ZLADIV
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, XERBLA, ZAXPY, ZDSCAL, ZTRSV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, DCONJG, DIMAG, MAX, MIN
      // ..
      // .. Statement Functions ..
      double             CABS1, CABS2;
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
      CABS2( ZDUM ) = ABS( DBLE( ZDUM ) / 2.D0 ) + ABS( DIMAG( ZDUM ) / 2.D0 )
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
         xerbla('ZLATRS', -INFO );
         RETURN
      }

      // Quick return if possible

      SCALE = ONE
      IF( N.EQ.0 ) RETURN

      // Determine machine dependent parameters to control overflow.

      SMLNUM = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      BIGNUM = ONE / SMLNUM

      if ( LSAME( NORMIN, 'N' ) ) {

         // Compute the 1-norm of each column, not including the diagonal.

         if ( UPPER ) {

            // A is upper triangular.

            DO 10 J = 1, N
               CNORM( J ) = DZASUM( J-1, A( 1, J ), 1 )
   10       CONTINUE
         } else {

            // A is lower triangular.

            DO 20 J = 1, N - 1
               CNORM( J ) = DZASUM( N-J, A( J+1, J ), 1 )
   20       CONTINUE
            CNORM( N ) = ZERO
         }
      }

      // Scale the column norms by TSCAL if the maximum element in CNORM is
      // greater than BIGNUM/2.

      IMAX = IDAMAX( N, CNORM, 1 )
      TMAX = CNORM( IMAX )
      if ( TMAX.LE.BIGNUM*HALF ) {
         TSCAL = ONE
      } else {

         // Avoid NaN generation if entries in CNORM exceed the
         // overflow threshold

         if ( TMAX.LE.DLAMCH('Overflow') ) {
            // Case 1: All entries in CNORM are valid floating-point numbers
            TSCAL = HALF / ( SMLNUM*TMAX )
            dscal(N, TSCAL, CNORM, 1 );
         } else {
            // Case 2: At least one column norm of A cannot be
            // represented as a floating-point number. Find the
            // maximum offdiagonal absolute value
            // max( |Re(A(I,J))|, |Im(A(I,J)| ). If this entry is
            // not +/- Infinity, use this value as TSCAL.
            TMAX = ZERO
            if ( UPPER ) {

               // A is upper triangular.

               DO J = 2, N
                  DO I = 1, J - 1
                     TMAX = MAX( TMAX, ABS( DBLE( A( I, J ) ) ), ABS( DIMAG(A ( I, J ) ) ) )
                  END DO
               END DO
            } else {

               // A is lower triangular.

               DO J = 1, N - 1
                  DO I = J + 1, N
                     TMAX = MAX( TMAX, ABS( DBLE( A( I, J ) ) ), ABS( DIMAG(A ( I, J ) ) ) )
                  END DO
               END DO
            }

            if ( TMAX.LE.DLAMCH('Overflow') ) {
               TSCAL = ONE / ( SMLNUM*TMAX )
               DO J = 1, N
                  if ( CNORM( J ).LE.DLAMCH('Overflow') ) {
                     CNORM( J ) = CNORM( J )*TSCAL
                  } else {
                     // Recompute the 1-norm of each column without
                     // introducing Infinity in the summation.
                     TSCAL = TWO * TSCAL
                     CNORM( J ) = ZERO
                     if ( UPPER ) {
                        DO I = 1, J - 1
                           CNORM( J ) = CNORM( J ) + TSCAL * CABS2( A( I, J ) )
                        END DO
                     } else {
                        DO I = J + 1, N
                           CNORM( J ) = CNORM( J ) + TSCAL * CABS2( A( I, J ) )
                        END DO
                     }
                     TSCAL = TSCAL * HALF
                  }
               END DO
            } else {
               // At least one entry of A is not a valid floating-point
               // entry. Rely on TRSV to propagate Inf and NaN.
               ztrsv(UPLO, TRANS, DIAG, N, A, LDA, X, 1 );
               RETURN
            }
         }
      }

      // Compute a bound on the computed solution vector to see if the
      // Level 2 BLAS routine ZTRSV can be used.

      XMAX = ZERO
      DO 30 J = 1, N
         XMAX = MAX( XMAX, CABS2( X( J ) ) )
   30 CONTINUE
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

               IF( GROW.LE.SMLNUM ) GO TO 60

               TJJS = A( J, J )
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
   40       CONTINUE
            GROW = XBND
         } else {

            // A is unit triangular.

            // Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.

            GROW = MIN( ONE, HALF / MAX( XBND, SMLNUM ) )
            DO 50 J = JFIRST, JLAST, JINC

               // Exit the loop if the growth factor is too small.

               IF( GROW.LE.SMLNUM ) GO TO 60

               // G(j) = G(j-1)*( 1 + CNORM(j) )

               GROW = GROW*( ONE / ( ONE+CNORM( J ) ) )
   50       CONTINUE
         }
   60    CONTINUE

      } else {

         // Compute the growth in A**T * x = b  or  A**H * x = b.

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

               IF( GROW.LE.SMLNUM ) GO TO 90

               // G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )

               XJ = ONE + CNORM( J )
               GROW = MIN( GROW, XBND / XJ )

               TJJS = A( J, J )
               TJJ = CABS1( TJJS )

               if ( TJJ.GE.SMLNUM ) {

                  // M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))

                  IF( XJ.GT.TJJ ) XBND = XBND*( TJJ / XJ )
               } else {

                  // M(j) could overflow, set XBND to 0.

                  XBND = ZERO
               }
   70       CONTINUE
            GROW = MIN( GROW, XBND )
         } else {

            // A is unit triangular.

            // Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.

            GROW = MIN( ONE, HALF / MAX( XBND, SMLNUM ) )
            DO 80 J = JFIRST, JLAST, JINC

               // Exit the loop if the growth factor is too small.

               IF( GROW.LE.SMLNUM ) GO TO 90

               // G(j) = ( 1 + CNORM(j) )*G(j-1)

               XJ = ONE + CNORM( J )
               GROW = GROW / XJ
   80       CONTINUE
         }
   90    CONTINUE
      }

      if ( ( GROW*TSCAL ).GT.SMLNUM ) {

         // Use the Level 2 BLAS solve if the reciprocal of the bound on
         // elements of X is not too small.

         ztrsv(UPLO, TRANS, DIAG, N, A, LDA, X, 1 );
      } else {

         // Use a Level 1 BLAS solve, scaling intermediate results.

         if ( XMAX.GT.BIGNUM*HALF ) {

            // Scale X so that its components are less than or equal to
            // BIGNUM in absolute value.

            SCALE = ( BIGNUM*HALF ) / XMAX
            zdscal(N, SCALE, X, 1 );
            XMAX = BIGNUM
         } else {
            XMAX = XMAX*TWO
         }

         if ( NOTRAN ) {

            // Solve A * x = b

            DO 120 J = JFIRST, JLAST, JINC

               // Compute x(j) = b(j) / A(j,j), scaling x if necessary.

               XJ = CABS1( X( J ) )
               if ( NOUNIT ) {
                  TJJS = A( J, J )*TSCAL
               } else {
                  TJJS = TSCAL
                  IF( TSCAL.EQ.ONE ) GO TO 110
               }
               TJJ = CABS1( TJJS )
               if ( TJJ.GT.SMLNUM ) {

                     // abs(A(j,j)) > SMLNUM:

                  if ( TJJ.LT.ONE ) {
                     if ( XJ.GT.TJJ*BIGNUM ) {

                           // Scale x by 1/b(j).

                        REC = ONE / XJ
                        zdscal(N, REC, X, 1 );
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     }
                  }
                  X( J ) = ZLADIV( X( J ), TJJS )
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
                     zdscal(N, REC, X, 1 );
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  }
                  X( J ) = ZLADIV( X( J ), TJJS )
                  XJ = CABS1( X( J ) )
               } else {

                     // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                     // scale = 0, and compute a solution to A*x = 0.

                  DO 100 I = 1, N
                     X( I ) = ZERO
  100             CONTINUE
                  X( J ) = ONE
                  XJ = ONE
                  SCALE = ZERO
                  XMAX = ZERO
               }
  110          CONTINUE

               // Scale x if necessary to avoid overflow when adding a
               // multiple of column j of A.

               if ( XJ.GT.ONE ) {
                  REC = ONE / XJ
                  if ( CNORM( J ).GT.( BIGNUM-XMAX )*REC ) {

                     // Scale x by 1/(2*abs(x(j))).

                     REC = REC*HALF
                     zdscal(N, REC, X, 1 );
                     SCALE = SCALE*REC
                  }
               } else if ( XJ*CNORM( J ).GT.( BIGNUM-XMAX ) ) {

                  // Scale x by 1/2.

                  zdscal(N, HALF, X, 1 );
                  SCALE = SCALE*HALF
               }

               if ( UPPER ) {
                  if ( J.GT.1 ) {

                     // Compute the update
                        // x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)

                     zaxpy(J-1, -X( J )*TSCAL, A( 1, J ), 1, X, 1 );
                     I = IZAMAX( J-1, X, 1 )
                     XMAX = CABS1( X( I ) )
                  }
               } else {
                  if ( J.LT.N ) {

                     // Compute the update
                        // x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)

                     zaxpy(N-J, -X( J )*TSCAL, A( J+1, J ), 1, X( J+1 ), 1 );
                     I = J + IZAMAX( N-J, X( J+1 ), 1 )
                     XMAX = CABS1( X( I ) )
                  }
               }
  120       CONTINUE

         } else if ( LSAME( TRANS, 'T' ) ) {

            // Solve A**T * x = b

            DO 170 J = JFIRST, JLAST, JINC

               // Compute x(j) = b(j) - sum A(k,j)*x(k).
                                     // k<>j

               XJ = CABS1( X( J ) )
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
                  TJJ = CABS1( TJJS )
                  if ( TJJ.GT.ONE ) {

                        // Divide by A(j,j) when scaling x if A(j,j) > 1.

                     REC = MIN( ONE, REC*TJJ )
                     USCAL = ZLADIV( USCAL, TJJS )
                  }
                  if ( REC.LT.ONE ) {
                     zdscal(N, REC, X, 1 );
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  }
               }

               CSUMJ = ZERO
               if ( USCAL.EQ.DCMPLX( ONE ) ) {

                  // If the scaling needed for A in the dot product is 1,
                  // call ZDOTU to perform the dot product.

                  if ( UPPER ) {
                     CSUMJ = ZDOTU( J-1, A( 1, J ), 1, X, 1 )
                  } else if ( J.LT.N ) {
                     CSUMJ = ZDOTU( N-J, A( J+1, J ), 1, X( J+1 ), 1 )
                  }
               } else {

                  // Otherwise, use in-line code for the dot product.

                  if ( UPPER ) {
                     DO 130 I = 1, J - 1
                        CSUMJ = CSUMJ + ( A( I, J )*USCAL )*X( I )
  130                CONTINUE
                  } else if ( J.LT.N ) {
                     DO 140 I = J + 1, N
                        CSUMJ = CSUMJ + ( A( I, J )*USCAL )*X( I )
  140                CONTINUE
                  }
               }

               if ( USCAL.EQ.DCMPLX( TSCAL ) ) {

                  // Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
                  // was not used to scale the dotproduct.

                  X( J ) = X( J ) - CSUMJ
                  XJ = CABS1( X( J ) )
                  if ( NOUNIT ) {
                     TJJS = A( J, J )*TSCAL
                  } else {
                     TJJS = TSCAL
                     IF( TSCAL.EQ.ONE ) GO TO 160
                  }

                     // Compute x(j) = x(j) / A(j,j), scaling if necessary.

                  TJJ = CABS1( TJJS )
                  if ( TJJ.GT.SMLNUM ) {

                        // abs(A(j,j)) > SMLNUM:

                     if ( TJJ.LT.ONE ) {
                        if ( XJ.GT.TJJ*BIGNUM ) {

                              // Scale X by 1/abs(x(j)).

                           REC = ONE / XJ
                           zdscal(N, REC, X, 1 );
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        }
                     }
                     X( J ) = ZLADIV( X( J ), TJJS )
                  } else if ( TJJ.GT.ZERO ) {

                        // 0 < abs(A(j,j)) <= SMLNUM:

                     if ( XJ.GT.TJJ*BIGNUM ) {

                           // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.

                        REC = ( TJJ*BIGNUM ) / XJ
                        zdscal(N, REC, X, 1 );
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     }
                     X( J ) = ZLADIV( X( J ), TJJS )
                  } else {

                        // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                        // scale = 0 and compute a solution to A**T *x = 0.

                     DO 150 I = 1, N
                        X( I ) = ZERO
  150                CONTINUE
                     X( J ) = ONE
                     SCALE = ZERO
                     XMAX = ZERO
                  }
  160             CONTINUE
               } else {

                  // Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
                  // product has already been divided by 1/A(j,j).

                  X( J ) = ZLADIV( X( J ), TJJS ) - CSUMJ
               }
               XMAX = MAX( XMAX, CABS1( X( J ) ) )
  170       CONTINUE

         } else {

            // Solve A**H * x = b

            DO 220 J = JFIRST, JLAST, JINC

               // Compute x(j) = b(j) - sum A(k,j)*x(k).
                                     // k<>j

               XJ = CABS1( X( J ) )
               USCAL = TSCAL
               REC = ONE / MAX( XMAX, ONE )
               if ( CNORM( J ).GT.( BIGNUM-XJ )*REC ) {

                  // If x(j) could overflow, scale x by 1/(2*XMAX).

                  REC = REC*HALF
                  if ( NOUNIT ) {
                     TJJS = DCONJG( A( J, J ) )*TSCAL
                  } else {
                     TJJS = TSCAL
                  }
                  TJJ = CABS1( TJJS )
                  if ( TJJ.GT.ONE ) {

                        // Divide by A(j,j) when scaling x if A(j,j) > 1.

                     REC = MIN( ONE, REC*TJJ )
                     USCAL = ZLADIV( USCAL, TJJS )
                  }
                  if ( REC.LT.ONE ) {
                     zdscal(N, REC, X, 1 );
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  }
               }

               CSUMJ = ZERO
               if ( USCAL.EQ.DCMPLX( ONE ) ) {

                  // If the scaling needed for A in the dot product is 1,
                  // call ZDOTC to perform the dot product.

                  if ( UPPER ) {
                     CSUMJ = ZDOTC( J-1, A( 1, J ), 1, X, 1 )
                  } else if ( J.LT.N ) {
                     CSUMJ = ZDOTC( N-J, A( J+1, J ), 1, X( J+1 ), 1 )
                  }
               } else {

                  // Otherwise, use in-line code for the dot product.

                  if ( UPPER ) {
                     DO 180 I = 1, J - 1
                        CSUMJ = CSUMJ + ( DCONJG( A( I, J ) )*USCAL )* X( I )
  180                CONTINUE
                  } else if ( J.LT.N ) {
                     DO 190 I = J + 1, N
                        CSUMJ = CSUMJ + ( DCONJG( A( I, J ) )*USCAL )* X( I )
  190                CONTINUE
                  }
               }

               if ( USCAL.EQ.DCMPLX( TSCAL ) ) {

                  // Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
                  // was not used to scale the dotproduct.

                  X( J ) = X( J ) - CSUMJ
                  XJ = CABS1( X( J ) )
                  if ( NOUNIT ) {
                     TJJS = DCONJG( A( J, J ) )*TSCAL
                  } else {
                     TJJS = TSCAL
                     IF( TSCAL.EQ.ONE ) GO TO 210
                  }

                     // Compute x(j) = x(j) / A(j,j), scaling if necessary.

                  TJJ = CABS1( TJJS )
                  if ( TJJ.GT.SMLNUM ) {

                        // abs(A(j,j)) > SMLNUM:

                     if ( TJJ.LT.ONE ) {
                        if ( XJ.GT.TJJ*BIGNUM ) {

                              // Scale X by 1/abs(x(j)).

                           REC = ONE / XJ
                           zdscal(N, REC, X, 1 );
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        }
                     }
                     X( J ) = ZLADIV( X( J ), TJJS )
                  } else if ( TJJ.GT.ZERO ) {

                        // 0 < abs(A(j,j)) <= SMLNUM:

                     if ( XJ.GT.TJJ*BIGNUM ) {

                           // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.

                        REC = ( TJJ*BIGNUM ) / XJ
                        zdscal(N, REC, X, 1 );
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     }
                     X( J ) = ZLADIV( X( J ), TJJS )
                  } else {

                        // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                        // scale = 0 and compute a solution to A**H *x = 0.

                     DO 200 I = 1, N
                        X( I ) = ZERO
  200                CONTINUE
                     X( J ) = ONE
                     SCALE = ZERO
                     XMAX = ZERO
                  }
  210             CONTINUE
               } else {

                  // Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
                  // product has already been divided by 1/A(j,j).

                  X( J ) = ZLADIV( X( J ), TJJS ) - CSUMJ
               }
               XMAX = MAX( XMAX, CABS1( X( J ) ) )
  220       CONTINUE
         }
         SCALE = SCALE / TSCAL
      }

      // Scale the column norms by 1/TSCAL for return.

      if ( TSCAL.NE.ONE ) {
         dscal(N, ONE / TSCAL, CNORM, 1 );
      }

      RETURN

      // End of ZLATRS

      }
