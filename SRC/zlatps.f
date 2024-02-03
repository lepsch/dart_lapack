      SUBROUTINE ZLATPS( UPLO, TRANS, DIAG, NORMIN, N, AP, X, SCALE, CNORM, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, NORMIN, TRANS, UPLO;
      int                INFO, N;
      double             SCALE;
      // ..
      // .. Array Arguments ..
      double             CNORM( * );
      COMPLEX*16         AP( * ), X( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, HALF, ONE, TWO;
      PARAMETER          ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
      // ..
      // .. Local Scalars ..
      bool               NOTRAN, NOUNIT, UPPER;
      int                I, IMAX, IP, J, JFIRST, JINC, JLAST, JLEN;
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
      // EXTERNAL DSCAL, XERBLA, ZAXPY, ZDSCAL, ZTPSV
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

      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( .NOT.LSAME( NORMIN, 'Y' ) .AND. .NOT. LSAME( NORMIN, 'N' ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZLATPS', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Determine machine dependent parameters to control overflow.

      SMLNUM = DLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SMLNUM
      SMLNUM = SMLNUM / DLAMCH( 'Precision' )
      BIGNUM = ONE / SMLNUM
      SCALE = ONE

      IF( LSAME( NORMIN, 'N' ) ) THEN

         // Compute the 1-norm of each column, not including the diagonal.

         IF( UPPER ) THEN

            // A is upper triangular.

            IP = 1
            DO 10 J = 1, N
               CNORM( J ) = DZASUM( J-1, AP( IP ), 1 )
               IP = IP + J
   10       CONTINUE
         ELSE

            // A is lower triangular.

            IP = 1
            DO 20 J = 1, N - 1
               CNORM( J ) = DZASUM( N-J, AP( IP+1 ), 1 )
               IP = IP + N - J + 1
   20       CONTINUE
            CNORM( N ) = ZERO
         END IF
      END IF

      // Scale the column norms by TSCAL if the maximum element in CNORM is
      // greater than BIGNUM/2.

      IMAX = IDAMAX( N, CNORM, 1 )
      TMAX = CNORM( IMAX )
      IF( TMAX.LE.BIGNUM*HALF ) THEN
         TSCAL = ONE
      ELSE
         TSCAL = HALF / ( SMLNUM*TMAX )
         CALL DSCAL( N, TSCAL, CNORM, 1 )
      END IF

      // Compute a bound on the computed solution vector to see if the
      // Level 2 BLAS routine ZTPSV can be used.

      XMAX = ZERO
      DO 30 J = 1, N
         XMAX = MAX( XMAX, CABS2( X( J ) ) )
   30 CONTINUE
      XBND = XMAX
      IF( NOTRAN ) THEN

         // Compute the growth in A * x = b.

         IF( UPPER ) THEN
            JFIRST = N
            JLAST = 1
            JINC = -1
         ELSE
            JFIRST = 1
            JLAST = N
            JINC = 1
         END IF

         IF( TSCAL.NE.ONE ) THEN
            GROW = ZERO
            GO TO 60
         END IF

         IF( NOUNIT ) THEN

            // A is non-unit triangular.

            // Compute GROW = 1/G(j) and XBND = 1/M(j).
            // Initially, G(0) = max{x(i), i=1,...,n}.

            GROW = HALF / MAX( XBND, SMLNUM )
            XBND = GROW
            IP = JFIRST*( JFIRST+1 ) / 2
            JLEN = N
            DO 40 J = JFIRST, JLAST, JINC

               // Exit the loop if the growth factor is too small.

               IF( GROW.LE.SMLNUM ) GO TO 60

               TJJS = AP( IP )
               TJJ = CABS1( TJJS )

               IF( TJJ.GE.SMLNUM ) THEN

                  // M(j) = G(j-1) / abs(A(j,j))

                  XBND = MIN( XBND, MIN( ONE, TJJ )*GROW )
               ELSE

                  // M(j) could overflow, set XBND to 0.

                  XBND = ZERO
               END IF

               IF( TJJ+CNORM( J ).GE.SMLNUM ) THEN

                  // G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )

                  GROW = GROW*( TJJ / ( TJJ+CNORM( J ) ) )
               ELSE

                  // G(j) could overflow, set GROW to 0.

                  GROW = ZERO
               END IF
               IP = IP + JINC*JLEN
               JLEN = JLEN - 1
   40       CONTINUE
            GROW = XBND
         ELSE

            // A is unit triangular.

            // Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.

            GROW = MIN( ONE, HALF / MAX( XBND, SMLNUM ) )
            DO 50 J = JFIRST, JLAST, JINC

               // Exit the loop if the growth factor is too small.

               IF( GROW.LE.SMLNUM ) GO TO 60

               // G(j) = G(j-1)*( 1 + CNORM(j) )

               GROW = GROW*( ONE / ( ONE+CNORM( J ) ) )
   50       CONTINUE
         END IF
   60    CONTINUE

      ELSE

         // Compute the growth in A**T * x = b  or  A**H * x = b.

         IF( UPPER ) THEN
            JFIRST = 1
            JLAST = N
            JINC = 1
         ELSE
            JFIRST = N
            JLAST = 1
            JINC = -1
         END IF

         IF( TSCAL.NE.ONE ) THEN
            GROW = ZERO
            GO TO 90
         END IF

         IF( NOUNIT ) THEN

            // A is non-unit triangular.

            // Compute GROW = 1/G(j) and XBND = 1/M(j).
            // Initially, M(0) = max{x(i), i=1,...,n}.

            GROW = HALF / MAX( XBND, SMLNUM )
            XBND = GROW
            IP = JFIRST*( JFIRST+1 ) / 2
            JLEN = 1
            DO 70 J = JFIRST, JLAST, JINC

               // Exit the loop if the growth factor is too small.

               IF( GROW.LE.SMLNUM ) GO TO 90

               // G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )

               XJ = ONE + CNORM( J )
               GROW = MIN( GROW, XBND / XJ )

               TJJS = AP( IP )
               TJJ = CABS1( TJJS )

               IF( TJJ.GE.SMLNUM ) THEN

                  // M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))

                  IF( XJ.GT.TJJ ) XBND = XBND*( TJJ / XJ )
               ELSE

                  // M(j) could overflow, set XBND to 0.

                  XBND = ZERO
               END IF
               JLEN = JLEN + 1
               IP = IP + JINC*JLEN
   70       CONTINUE
            GROW = MIN( GROW, XBND )
         ELSE

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
         END IF
   90    CONTINUE
      END IF

      IF( ( GROW*TSCAL ).GT.SMLNUM ) THEN

         // Use the Level 2 BLAS solve if the reciprocal of the bound on
         // elements of X is not too small.

         CALL ZTPSV( UPLO, TRANS, DIAG, N, AP, X, 1 )
      ELSE

         // Use a Level 1 BLAS solve, scaling intermediate results.

         IF( XMAX.GT.BIGNUM*HALF ) THEN

            // Scale X so that its components are less than or equal to
            // BIGNUM in absolute value.

            SCALE = ( BIGNUM*HALF ) / XMAX
            CALL ZDSCAL( N, SCALE, X, 1 )
            XMAX = BIGNUM
         ELSE
            XMAX = XMAX*TWO
         END IF

         IF( NOTRAN ) THEN

            // Solve A * x = b

            IP = JFIRST*( JFIRST+1 ) / 2
            DO 120 J = JFIRST, JLAST, JINC

               // Compute x(j) = b(j) / A(j,j), scaling x if necessary.

               XJ = CABS1( X( J ) )
               IF( NOUNIT ) THEN
                  TJJS = AP( IP )*TSCAL
               ELSE
                  TJJS = TSCAL
                  IF( TSCAL.EQ.ONE ) GO TO 110
               END IF
               TJJ = CABS1( TJJS )
               IF( TJJ.GT.SMLNUM ) THEN

                     // abs(A(j,j)) > SMLNUM:

                  IF( TJJ.LT.ONE ) THEN
                     IF( XJ.GT.TJJ*BIGNUM ) THEN

                           // Scale x by 1/b(j).

                        REC = ONE / XJ
                        CALL ZDSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
                  X( J ) = ZLADIV( X( J ), TJJS )
                  XJ = CABS1( X( J ) )
               ELSE IF( TJJ.GT.ZERO ) THEN

                     // 0 < abs(A(j,j)) <= SMLNUM:

                  IF( XJ.GT.TJJ*BIGNUM ) THEN

                        // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
                       t // o avoid overflow when dividing by A(j,j).

                     REC = ( TJJ*BIGNUM ) / XJ
                     IF( CNORM( J ).GT.ONE ) THEN

                           // Scale by 1/CNORM(j) to avoid overflow when
                           // multiplying x(j) times column j.

                        REC = REC / CNORM( J )
                     END IF
                     CALL ZDSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
                  X( J ) = ZLADIV( X( J ), TJJS )
                  XJ = CABS1( X( J ) )
               ELSE

                     // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                     // scale = 0, and compute a solution to A*x = 0.

                  DO 100 I = 1, N
                     X( I ) = ZERO
  100             CONTINUE
                  X( J ) = ONE
                  XJ = ONE
                  SCALE = ZERO
                  XMAX = ZERO
               END IF
  110          CONTINUE

               // Scale x if necessary to avoid overflow when adding a
               // multiple of column j of A.

               IF( XJ.GT.ONE ) THEN
                  REC = ONE / XJ
                  IF( CNORM( J ).GT.( BIGNUM-XMAX )*REC ) THEN

                     // Scale x by 1/(2*abs(x(j))).

                     REC = REC*HALF
                     CALL ZDSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                  END IF
               ELSE IF( XJ*CNORM( J ).GT.( BIGNUM-XMAX ) ) THEN

                  // Scale x by 1/2.

                  CALL ZDSCAL( N, HALF, X, 1 )
                  SCALE = SCALE*HALF
               END IF

               IF( UPPER ) THEN
                  IF( J.GT.1 ) THEN

                     // Compute the update
                        // x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)

                     CALL ZAXPY( J-1, -X( J )*TSCAL, AP( IP-J+1 ), 1, X, 1 )
                     I = IZAMAX( J-1, X, 1 )
                     XMAX = CABS1( X( I ) )
                  END IF
                  IP = IP - J
               ELSE
                  IF( J.LT.N ) THEN

                     // Compute the update
                        // x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)

                     CALL ZAXPY( N-J, -X( J )*TSCAL, AP( IP+1 ), 1, X( J+1 ), 1 )
                     I = J + IZAMAX( N-J, X( J+1 ), 1 )
                     XMAX = CABS1( X( I ) )
                  END IF
                  IP = IP + N - J + 1
               END IF
  120       CONTINUE

         ELSE IF( LSAME( TRANS, 'T' ) ) THEN

            // Solve A**T * x = b

            IP = JFIRST*( JFIRST+1 ) / 2
            JLEN = 1
            DO 170 J = JFIRST, JLAST, JINC

               // Compute x(j) = b(j) - sum A(k,j)*x(k).
                                     // k<>j

               XJ = CABS1( X( J ) )
               USCAL = TSCAL
               REC = ONE / MAX( XMAX, ONE )
               IF( CNORM( J ).GT.( BIGNUM-XJ )*REC ) THEN

                  // If x(j) could overflow, scale x by 1/(2*XMAX).

                  REC = REC*HALF
                  IF( NOUNIT ) THEN
                     TJJS = AP( IP )*TSCAL
                  ELSE
                     TJJS = TSCAL
                  END IF
                  TJJ = CABS1( TJJS )
                  IF( TJJ.GT.ONE ) THEN

                        // Divide by A(j,j) when scaling x if A(j,j) > 1.

                     REC = MIN( ONE, REC*TJJ )
                     USCAL = ZLADIV( USCAL, TJJS )
                  END IF
                  IF( REC.LT.ONE ) THEN
                     CALL ZDSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
               END IF

               CSUMJ = ZERO
               IF( USCAL.EQ.DCMPLX( ONE ) ) THEN

                  // If the scaling needed for A in the dot product is 1,
                  // call ZDOTU to perform the dot product.

                  IF( UPPER ) THEN
                     CSUMJ = ZDOTU( J-1, AP( IP-J+1 ), 1, X, 1 )
                  ELSE IF( J.LT.N ) THEN
                     CSUMJ = ZDOTU( N-J, AP( IP+1 ), 1, X( J+1 ), 1 )
                  END IF
               ELSE

                  // Otherwise, use in-line code for the dot product.

                  IF( UPPER ) THEN
                     DO 130 I = 1, J - 1
                        CSUMJ = CSUMJ + ( AP( IP-J+I )*USCAL )*X( I )
  130                CONTINUE
                  ELSE IF( J.LT.N ) THEN
                     DO 140 I = 1, N - J
                        CSUMJ = CSUMJ + ( AP( IP+I )*USCAL )*X( J+I )
  140                CONTINUE
                  END IF
               END IF

               IF( USCAL.EQ.DCMPLX( TSCAL ) ) THEN

                  // Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
                  // was not used to scale the dotproduct.

                  X( J ) = X( J ) - CSUMJ
                  XJ = CABS1( X( J ) )
                  IF( NOUNIT ) THEN

                     // Compute x(j) = x(j) / A(j,j), scaling if necessary.

                     TJJS = AP( IP )*TSCAL
                  ELSE
                     TJJS = TSCAL
                     IF( TSCAL.EQ.ONE ) GO TO 160
                  END IF
                  TJJ = CABS1( TJJS )
                  IF( TJJ.GT.SMLNUM ) THEN

                        // abs(A(j,j)) > SMLNUM:

                     IF( TJJ.LT.ONE ) THEN
                        IF( XJ.GT.TJJ*BIGNUM ) THEN

                              // Scale X by 1/abs(x(j)).

                           REC = ONE / XJ
                           CALL ZDSCAL( N, REC, X, 1 )
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        END IF
                     END IF
                     X( J ) = ZLADIV( X( J ), TJJS )
                  ELSE IF( TJJ.GT.ZERO ) THEN

                        // 0 < abs(A(j,j)) <= SMLNUM:

                     IF( XJ.GT.TJJ*BIGNUM ) THEN

                           // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.

                        REC = ( TJJ*BIGNUM ) / XJ
                        CALL ZDSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                     X( J ) = ZLADIV( X( J ), TJJS )
                  ELSE

                        // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                        // scale = 0 and compute a solution to A**T *x = 0.

                     DO 150 I = 1, N
                        X( I ) = ZERO
  150                CONTINUE
                     X( J ) = ONE
                     SCALE = ZERO
                     XMAX = ZERO
                  END IF
  160             CONTINUE
               ELSE

                  // Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
                  // product has already been divided by 1/A(j,j).

                  X( J ) = ZLADIV( X( J ), TJJS ) - CSUMJ
               END IF
               XMAX = MAX( XMAX, CABS1( X( J ) ) )
               JLEN = JLEN + 1
               IP = IP + JINC*JLEN
  170       CONTINUE

         ELSE

            // Solve A**H * x = b

            IP = JFIRST*( JFIRST+1 ) / 2
            JLEN = 1
            DO 220 J = JFIRST, JLAST, JINC

               // Compute x(j) = b(j) - sum A(k,j)*x(k).
                                     // k<>j

               XJ = CABS1( X( J ) )
               USCAL = TSCAL
               REC = ONE / MAX( XMAX, ONE )
               IF( CNORM( J ).GT.( BIGNUM-XJ )*REC ) THEN

                  // If x(j) could overflow, scale x by 1/(2*XMAX).

                  REC = REC*HALF
                  IF( NOUNIT ) THEN
                     TJJS = DCONJG( AP( IP ) )*TSCAL
                  ELSE
                     TJJS = TSCAL
                  END IF
                  TJJ = CABS1( TJJS )
                  IF( TJJ.GT.ONE ) THEN

                        // Divide by A(j,j) when scaling x if A(j,j) > 1.

                     REC = MIN( ONE, REC*TJJ )
                     USCAL = ZLADIV( USCAL, TJJS )
                  END IF
                  IF( REC.LT.ONE ) THEN
                     CALL ZDSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
               END IF

               CSUMJ = ZERO
               IF( USCAL.EQ.DCMPLX( ONE ) ) THEN

                  // If the scaling needed for A in the dot product is 1,
                  // call ZDOTC to perform the dot product.

                  IF( UPPER ) THEN
                     CSUMJ = ZDOTC( J-1, AP( IP-J+1 ), 1, X, 1 )
                  ELSE IF( J.LT.N ) THEN
                     CSUMJ = ZDOTC( N-J, AP( IP+1 ), 1, X( J+1 ), 1 )
                  END IF
               ELSE

                  // Otherwise, use in-line code for the dot product.

                  IF( UPPER ) THEN
                     DO 180 I = 1, J - 1
                        CSUMJ = CSUMJ + ( DCONJG( AP( IP-J+I ) )*USCAL ) *X( I )
  180                CONTINUE
                  ELSE IF( J.LT.N ) THEN
                     DO 190 I = 1, N - J
                        CSUMJ = CSUMJ + ( DCONJG( AP( IP+I ) )*USCAL )* X( J+I )
  190                CONTINUE
                  END IF
               END IF

               IF( USCAL.EQ.DCMPLX( TSCAL ) ) THEN

                  // Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
                  // was not used to scale the dotproduct.

                  X( J ) = X( J ) - CSUMJ
                  XJ = CABS1( X( J ) )
                  IF( NOUNIT ) THEN

                     // Compute x(j) = x(j) / A(j,j), scaling if necessary.

                     TJJS = DCONJG( AP( IP ) )*TSCAL
                  ELSE
                     TJJS = TSCAL
                     IF( TSCAL.EQ.ONE ) GO TO 210
                  END IF
                  TJJ = CABS1( TJJS )
                  IF( TJJ.GT.SMLNUM ) THEN

                        // abs(A(j,j)) > SMLNUM:

                     IF( TJJ.LT.ONE ) THEN
                        IF( XJ.GT.TJJ*BIGNUM ) THEN

                              // Scale X by 1/abs(x(j)).

                           REC = ONE / XJ
                           CALL ZDSCAL( N, REC, X, 1 )
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        END IF
                     END IF
                     X( J ) = ZLADIV( X( J ), TJJS )
                  ELSE IF( TJJ.GT.ZERO ) THEN

                        // 0 < abs(A(j,j)) <= SMLNUM:

                     IF( XJ.GT.TJJ*BIGNUM ) THEN

                           // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.

                        REC = ( TJJ*BIGNUM ) / XJ
                        CALL ZDSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                     X( J ) = ZLADIV( X( J ), TJJS )
                  ELSE

                        // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                        // scale = 0 and compute a solution to A**H *x = 0.

                     DO 200 I = 1, N
                        X( I ) = ZERO
  200                CONTINUE
                     X( J ) = ONE
                     SCALE = ZERO
                     XMAX = ZERO
                  END IF
  210             CONTINUE
               ELSE

                  // Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
                  // product has already been divided by 1/A(j,j).

                  X( J ) = ZLADIV( X( J ), TJJS ) - CSUMJ
               END IF
               XMAX = MAX( XMAX, CABS1( X( J ) ) )
               JLEN = JLEN + 1
               IP = IP + JINC*JLEN
  220       CONTINUE
         END IF
         SCALE = SCALE / TSCAL
      END IF

      // Scale the column norms by 1/TSCAL for return.

      IF( TSCAL.NE.ONE ) THEN
         CALL DSCAL( N, ONE / TSCAL, CNORM, 1 )
      END IF

      RETURN

      // End of ZLATPS

      END
