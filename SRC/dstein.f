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
      const              ZERO = 0.0D+0, ONE = 1.0D+0, TEN = 1.0D+1, ODM3 = 1.0D-3, ODM1 = 1.0D-1 ;
      int                MAXITS, EXTRA;
      const              MAXITS = 5, EXTRA = 2 ;
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
      for (I = 1; I <= M; I++) { // 10
         IFAIL( I ) = 0
      } // 10

      if ( N.LT.0 ) {
         INFO = -1
      } else if ( M.LT.0 .OR. M.GT.N ) {
         INFO = -4
      } else if ( LDZ.LT.MAX( 1, N ) ) {
         INFO = -9
      } else {
         for (J = 2; J <= M; J++) { // 20
            if ( IBLOCK( J ).LT.IBLOCK( J-1 ) ) {
               INFO = -6
               GO TO 30
            }
            if ( IBLOCK( J ).EQ.IBLOCK( J-1 ) .AND. W( J ).LT.W( J-1 ) ) {
               INFO = -5
               GO TO 30
            }
         } // 20
         } // 30
      }

      if ( INFO.NE.0 ) {
         xerbla('DSTEIN', -INFO );
         RETURN
      }

      // Quick return if possible

      if ( N.EQ.0 .OR. M.EQ.0 ) {
         RETURN
      } else if ( N.EQ.1 ) {
         Z( 1, 1 ) = ONE
         RETURN
      }

      // Get machine constants.

      EPS = DLAMCH( 'Precision' )

      // Initialize seed for random number generator DLARNV.

      for (I = 1; I <= 4; I++) { // 40
         ISEED( I ) = 1
      } // 40

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

         if ( NBLK.EQ.1 ) {
            B1 = 1
         } else {
            B1 = ISPLIT( NBLK-1 ) + 1
         }
         BN = ISPLIT( NBLK )
         BLKSIZ = BN - B1 + 1
         IF( BLKSIZ.EQ.1 ) GO TO 60
         GPIND = J1

         // Compute reorthogonalization criterion and stopping criterion.

         ONENRM = ABS( D( B1 ) ) + ABS( E( B1 ) )
         ONENRM = MAX( ONENRM, ABS( D( BN ) )+ABS( E( BN-1 ) ) )
         DO 50 I = B1 + 1, BN - 1
            ONENRM = MAX( ONENRM, ABS( D( I ) )+ABS( E( I-1 ) )+ ABS( E( I ) ) )
         } // 50
         ORTOL = ODM3*ONENRM

         DTPCRT = SQRT( ODM1 / BLKSIZ )

         // Loop through eigenvalues of block nblk.

         } // 60
         JBLK = 0
         for (J = J1; J <= M; J++) { // 150
            if ( IBLOCK( J ).NE.NBLK ) {
               J1 = J
               GO TO 160
            }
            JBLK = JBLK + 1
            XJ = W( J )

            // Skip all the work if the block size is one.

            if ( BLKSIZ.EQ.1 ) {
               WORK( INDRV1+1 ) = ONE
               GO TO 120
            }

            // If eigenvalues j and j-1 are too close, add a relatively
            // small perturbation.

            if ( JBLK.GT.1 ) {
               EPS1 = ABS( EPS*XJ )
               PERTOL = TEN*EPS1
               SEP = XJ - XJM
               IF( SEP.LT.PERTOL ) XJ = XJM + PERTOL
            }

            ITS = 0
            NRMCHK = 0

            // Get random starting vector.

            dlarnv(2, ISEED, BLKSIZ, WORK( INDRV1+1 ) );

            // Copy the matrix T so it won't be destroyed in factorization.

            dcopy(BLKSIZ, D( B1 ), 1, WORK( INDRV4+1 ), 1 );
            dcopy(BLKSIZ-1, E( B1 ), 1, WORK( INDRV2+2 ), 1 );
            dcopy(BLKSIZ-1, E( B1 ), 1, WORK( INDRV3+1 ), 1 );

            // Compute LU factors with partial pivoting  ( PT = LU )

            TOL = ZERO
            dlagtf(BLKSIZ, WORK( INDRV4+1 ), XJ, WORK( INDRV2+2 ), WORK( INDRV3+1 ), TOL, WORK( INDRV5+1 ), IWORK, IINFO );

            // Update iteration count.

            } // 70
            ITS = ITS + 1
            IF( ITS.GT.MAXITS ) GO TO 100

            // Normalize and scale the righthand side vector Pb.

            JMAX = IDAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
            SCL = BLKSIZ*ONENRM*MAX( EPS, ABS( WORK( INDRV4+BLKSIZ ) ) ) / ABS( WORK( INDRV1+JMAX ) )
            dscal(BLKSIZ, SCL, WORK( INDRV1+1 ), 1 );

            // Solve the system LU = Pb.

            dlagts(-1, BLKSIZ, WORK( INDRV4+1 ), WORK( INDRV2+2 ), WORK( INDRV3+1 ), WORK( INDRV5+1 ), IWORK, WORK( INDRV1+1 ), TOL, IINFO );

            // Reorthogonalize by modified Gram-Schmidt if eigenvalues are
            // close enough.

            IF( JBLK.EQ.1 ) GO TO 90             IF( ABS( XJ-XJM ).GT.ORTOL ) GPIND = J
            if ( GPIND.NE.J ) {
               DO 80 I = GPIND, J - 1
                  ZTR = -DDOT( BLKSIZ, WORK( INDRV1+1 ), 1, Z( B1, I ), 1 )                   CALL DAXPY( BLKSIZ, ZTR, Z( B1, I ), 1, WORK( INDRV1+1 ), 1 )
               } // 80
            }

            // Check the infinity norm of the iterate.

            } // 90
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

            } // 100
            INFO = INFO + 1
            IFAIL( INFO ) = J

            // Accept iterate as jth eigenvector.

            } // 110
            SCL = ONE / DNRM2( BLKSIZ, WORK( INDRV1+1 ), 1 )
            JMAX = IDAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
            IF( WORK( INDRV1+JMAX ).LT.ZERO ) SCL = -SCL
            dscal(BLKSIZ, SCL, WORK( INDRV1+1 ), 1 );
            } // 120
            for (I = 1; I <= N; I++) { // 130
               Z( I, J ) = ZERO
            } // 130
            for (I = 1; I <= BLKSIZ; I++) { // 140
               Z( B1+I-1, J ) = WORK( INDRV1+I )
            } // 140

            // Save the shift to check eigenvalue spacing at next
            // iteration.

            XJM = XJ

         } // 150
      } // 160

      RETURN

      // End of DSTEIN

      }
