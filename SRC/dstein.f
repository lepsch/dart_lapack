      SUBROUTINE DSTEIN( N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK, IWORK, IFAIL, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDZ, M, N;
      // ..
      // .. Array Arguments ..
      int                IBLOCK( * ), IFAIL( * ), ISPLIT( * ), IWORK( * );
      double             D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TEN, ODM3, ODM1;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TEN = 1.0D+1, ODM3 = 1.0D-3, ODM1 = 1.0D-1 )
      int                MAXITS, EXTRA;
      PARAMETER          ( MAXITS = 5, EXTRA = 2 )
      // ..
      // .. Local Scalars ..
      int                B1, BLKSIZ, BN, GPIND, I, IINFO, INDRV1, INDRV2, INDRV3, INDRV4, INDRV5, ITS, J, J1, JBLK, JMAX, NBLK, NRMCHK;
      double             DTPCRT, EPS, EPS1, NRM, ONENRM, ORTOL, PERTOL, SCL, SEP, TOL, XJ, XJM, ZTR;
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 );
      // ..
      // .. External Functions ..
      int                IDAMAX;
      double             DDOT, DLAMCH, DNRM2;
      // EXTERNAL IDAMAX, DDOT, DLAMCH, DNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DLAGTF, DLAGTS, DLARNV, DSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      DO 10 I = 1, M
         IFAIL( I ) = 0
   10 CONTINUE

      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 .OR. M.GT.N ) THEN
         INFO = -4
      ELSE IF( LDZ.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE
         DO 20 J = 2, M
            IF( IBLOCK( J ).LT.IBLOCK( J-1 ) ) THEN
               INFO = -6
               GO TO 30
            END IF
            IF( IBLOCK( J ).EQ.IBLOCK( J-1 ) .AND. W( J ).LT.W( J-1 ) ) THEN
               INFO = -5
               GO TO 30
            END IF
   20    CONTINUE
   30    CONTINUE
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSTEIN', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 .OR. M.EQ.0 ) THEN
         RETURN
      ELSE IF( N.EQ.1 ) THEN
         Z( 1, 1 ) = ONE
         RETURN
      END IF

      // Get machine constants.

      EPS = DLAMCH( 'Precision' )

      // Initialize seed for random number generator DLARNV.

      DO 40 I = 1, 4
         ISEED( I ) = 1
   40 CONTINUE

      // Initialize pointers.

      INDRV1 = 0
      INDRV2 = INDRV1 + N
      INDRV3 = INDRV2 + N
      INDRV4 = INDRV3 + N
      INDRV5 = INDRV4 + N

      // Compute eigenvectors of matrix blocks.

      J1 = 1
      DO 160 NBLK = 1, IBLOCK( M )

         // Find starting and ending indices of block nblk.

         IF( NBLK.EQ.1 ) THEN
            B1 = 1
         ELSE
            B1 = ISPLIT( NBLK-1 ) + 1
         END IF
         BN = ISPLIT( NBLK )
         BLKSIZ = BN - B1 + 1
         IF( BLKSIZ.EQ.1 ) GO TO 60
         GPIND = J1

         // Compute reorthogonalization criterion and stopping criterion.

         ONENRM = ABS( D( B1 ) ) + ABS( E( B1 ) )
         ONENRM = MAX( ONENRM, ABS( D( BN ) )+ABS( E( BN-1 ) ) )
         DO 50 I = B1 + 1, BN - 1
            ONENRM = MAX( ONENRM, ABS( D( I ) )+ABS( E( I-1 ) )+ ABS( E( I ) ) )
   50    CONTINUE
         ORTOL = ODM3*ONENRM

         DTPCRT = SQRT( ODM1 / BLKSIZ )

         // Loop through eigenvalues of block nblk.

   60    CONTINUE
         JBLK = 0
         DO 150 J = J1, M
            IF( IBLOCK( J ).NE.NBLK ) THEN
               J1 = J
               GO TO 160
            END IF
            JBLK = JBLK + 1
            XJ = W( J )

            // Skip all the work if the block size is one.

            IF( BLKSIZ.EQ.1 ) THEN
               WORK( INDRV1+1 ) = ONE
               GO TO 120
            END IF

            // If eigenvalues j and j-1 are too close, add a relatively
            // small perturbation.

            IF( JBLK.GT.1 ) THEN
               EPS1 = ABS( EPS*XJ )
               PERTOL = TEN*EPS1
               SEP = XJ - XJM
               IF( SEP.LT.PERTOL ) XJ = XJM + PERTOL
            END IF

            ITS = 0
            NRMCHK = 0

            // Get random starting vector.

            CALL DLARNV( 2, ISEED, BLKSIZ, WORK( INDRV1+1 ) )

            // Copy the matrix T so it won't be destroyed in factorization.

            CALL DCOPY( BLKSIZ, D( B1 ), 1, WORK( INDRV4+1 ), 1 )
            CALL DCOPY( BLKSIZ-1, E( B1 ), 1, WORK( INDRV2+2 ), 1 )
            CALL DCOPY( BLKSIZ-1, E( B1 ), 1, WORK( INDRV3+1 ), 1 )

            // Compute LU factors with partial pivoting  ( PT = LU )

            TOL = ZERO
            CALL DLAGTF( BLKSIZ, WORK( INDRV4+1 ), XJ, WORK( INDRV2+2 ), WORK( INDRV3+1 ), TOL, WORK( INDRV5+1 ), IWORK, IINFO )

            // Update iteration count.

   70       CONTINUE
            ITS = ITS + 1
            IF( ITS.GT.MAXITS ) GO TO 100

            // Normalize and scale the righthand side vector Pb.

            JMAX = IDAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
            SCL = BLKSIZ*ONENRM*MAX( EPS, ABS( WORK( INDRV4+BLKSIZ ) ) ) / ABS( WORK( INDRV1+JMAX ) )
            CALL DSCAL( BLKSIZ, SCL, WORK( INDRV1+1 ), 1 )

            // Solve the system LU = Pb.

            CALL DLAGTS( -1, BLKSIZ, WORK( INDRV4+1 ), WORK( INDRV2+2 ), WORK( INDRV3+1 ), WORK( INDRV5+1 ), IWORK, WORK( INDRV1+1 ), TOL, IINFO )

            // Reorthogonalize by modified Gram-Schmidt if eigenvalues are
            // close enough.

            IF( JBLK.EQ.1 ) GO TO 90             IF( ABS( XJ-XJM ).GT.ORTOL ) GPIND = J
            IF( GPIND.NE.J ) THEN
               DO 80 I = GPIND, J - 1
                  ZTR = -DDOT( BLKSIZ, WORK( INDRV1+1 ), 1, Z( B1, I ), 1 )                   CALL DAXPY( BLKSIZ, ZTR, Z( B1, I ), 1, WORK( INDRV1+1 ), 1 )
   80          CONTINUE
            END IF

            // Check the infinity norm of the iterate.

   90       CONTINUE
            JMAX = IDAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
            NRM = ABS( WORK( INDRV1+JMAX ) )

            // Continue for additional iterations after norm reaches
            // stopping criterion.

            IF( NRM.LT.DTPCRT ) GO TO 70
            NRMCHK = NRMCHK + 1
            IF( NRMCHK.LT.EXTRA+1 ) GO TO 70

            GO TO 110

            // If stopping criterion was not satisfied, update info and
            // store eigenvector number in array ifail.

  100       CONTINUE
            INFO = INFO + 1
            IFAIL( INFO ) = J

            // Accept iterate as jth eigenvector.

  110       CONTINUE
            SCL = ONE / DNRM2( BLKSIZ, WORK( INDRV1+1 ), 1 )
            JMAX = IDAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
            IF( WORK( INDRV1+JMAX ).LT.ZERO ) SCL = -SCL
            CALL DSCAL( BLKSIZ, SCL, WORK( INDRV1+1 ), 1 )
  120       CONTINUE
            DO 130 I = 1, N
               Z( I, J ) = ZERO
  130       CONTINUE
            DO 140 I = 1, BLKSIZ
               Z( B1+I-1, J ) = WORK( INDRV1+I )
  140       CONTINUE

            // Save the shift to check eigenvalue spacing at next
            // iteration.

            XJM = XJ

  150    CONTINUE
  160 CONTINUE

      RETURN

      // End of DSTEIN

      END
