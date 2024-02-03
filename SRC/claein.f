      SUBROUTINE CLAEIN( RIGHTV, NOINIT, N, H, LDH, W, V, B, LDB, RWORK, EPS3, SMLNUM, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               NOINIT, RIGHTV;
      int                INFO, LDB, LDH, N;
      REAL               EPS3, SMLNUM
      COMPLEX            W
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            B( LDB, * ), H( LDH, * ), V( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, TENTH
      const              ONE = 1.0E+0, TENTH = 1.0E-1 ;
      COMPLEX            ZERO
      const              ZERO = ( 0.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      String             NORMIN, TRANS;
      int                I, IERR, ITS, J;
      REAL               GROWTO, NRMSML, ROOTN, RTEMP, SCALE, VNORM
      COMPLEX            CDUM, EI, EJ, TEMP, X
      // ..
      // .. External Functions ..
      int                ICAMAX;
      REAL               SCASUM, SCNRM2
      COMPLEX            CLADIV
      // EXTERNAL ICAMAX, SCASUM, SCNRM2, CLADIV
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLATRS, CSSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, REAL, SQRT
      // ..
      // .. Statement Functions ..
      REAL               CABS1
      // ..
      // .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
      // ..
      // .. Executable Statements ..

      INFO = 0

      // GROWTO is the threshold used in the acceptance test for an
      // eigenvector.

      ROOTN = SQRT( REAL( N ) )
      GROWTO = TENTH / ROOTN
      NRMSML = MAX( ONE, EPS3*ROOTN )*SMLNUM

      // Form B = H - W*I (except that the subdiagonal elements are not
      // stored).

      for (J = 1; J <= N; J++) { // 20
         DO 10 I = 1, J - 1
            B( I, J ) = H( I, J )
   10    CONTINUE
         B( J, J ) = H( J, J ) - W
   20 CONTINUE

      if ( NOINIT ) {

         // Initialize V.

         for (I = 1; I <= N; I++) { // 30
            V( I ) = EPS3
   30    CONTINUE
      } else {

         // Scale supplied initial vector.

         VNORM = SCNRM2( N, V, 1 )
         csscal(N, ( EPS3*ROOTN ) / MAX( VNORM, NRMSML ), V, 1 );
      }

      if ( RIGHTV ) {

         // LU decomposition with partial pivoting of B, replacing zero
         // pivots by EPS3.

         DO 60 I = 1, N - 1
            EI = H( I+1, I )
            if ( CABS1( B( I, I ) ).LT.CABS1( EI ) ) {

               // Interchange rows and eliminate.

               X = CLADIV( B( I, I ), EI )
               B( I, I ) = EI
               DO 40 J = I + 1, N
                  TEMP = B( I+1, J )
                  B( I+1, J ) = B( I, J ) - X*TEMP
                  B( I, J ) = TEMP
   40          CONTINUE
            } else {

               // Eliminate without interchange.

               IF( B( I, I ).EQ.ZERO ) B( I, I ) = EPS3
               X = CLADIV( EI, B( I, I ) )
               if ( X.NE.ZERO ) {
                  DO 50 J = I + 1, N
                     B( I+1, J ) = B( I+1, J ) - X*B( I, J )
   50             CONTINUE
               }
            }
   60    CONTINUE
         IF( B( N, N ).EQ.ZERO ) B( N, N ) = EPS3

         TRANS = 'N'

      } else {

         // UL decomposition with partial pivoting of B, replacing zero
         // pivots by EPS3.

         DO 90 J = N, 2, -1
            EJ = H( J, J-1 )
            if ( CABS1( B( J, J ) ).LT.CABS1( EJ ) ) {

               // Interchange columns and eliminate.

               X = CLADIV( B( J, J ), EJ )
               B( J, J ) = EJ
               DO 70 I = 1, J - 1
                  TEMP = B( I, J-1 )
                  B( I, J-1 ) = B( I, J ) - X*TEMP
                  B( I, J ) = TEMP
   70          CONTINUE
            } else {

               // Eliminate without interchange.

               IF( B( J, J ).EQ.ZERO ) B( J, J ) = EPS3
               X = CLADIV( EJ, B( J, J ) )
               if ( X.NE.ZERO ) {
                  DO 80 I = 1, J - 1
                     B( I, J-1 ) = B( I, J-1 ) - X*B( I, J )
   80             CONTINUE
               }
            }
   90    CONTINUE
         IF( B( 1, 1 ).EQ.ZERO ) B( 1, 1 ) = EPS3

         TRANS = 'C'

      }

      NORMIN = 'N'
      for (ITS = 1; ITS <= N; ITS++) { // 110

         // Solve U*x = scale*v for a right eigenvector
           // or U**H *x = scale*v for a left eigenvector,
         // overwriting x on v.

         clatrs('Upper', TRANS, 'Nonunit', NORMIN, N, B, LDB, V, SCALE, RWORK, IERR );
         NORMIN = 'Y'

         // Test for sufficient growth in the norm of v.

         VNORM = SCASUM( N, V, 1 )
         IF( VNORM.GE.GROWTO*SCALE ) GO TO 120

         // Choose new orthogonal starting vector and try again.

         RTEMP = EPS3 / ( ROOTN+ONE )
         V( 1 ) = EPS3
         for (I = 2; I <= N; I++) { // 100
            V( I ) = RTEMP
  100    CONTINUE
         V( N-ITS+1 ) = V( N-ITS+1 ) - EPS3*ROOTN
  110 CONTINUE

      // Failure to find eigenvector in N iterations.

      INFO = 1

  120 CONTINUE

      // Normalize eigenvector.

      I = ICAMAX( N, V, 1 )
      csscal(N, ONE / CABS1( V( I ) ), V, 1 );

      RETURN

      // End of CLAEIN

      }
