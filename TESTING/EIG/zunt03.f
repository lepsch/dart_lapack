      SUBROUTINE ZUNT03( RC, MU, MV, N, K, U, LDU, V, LDV, WORK, LWORK, RWORK, RESULT, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      List<String>       RC;
      int                INFO, K, LDU, LDV, LWORK, MU, MV, N;
      double             RESULT;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      COMPLEX*16         U( LDU, * ), V( LDV, * ), WORK( * )
      // ..

*  =====================================================================


      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      int                I, IRC, J, LMX;
      double             RES1, RES2, ULP;
      COMPLEX*16         S, SU, SV
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IZAMAX;
      double             DLAMCH;
      // EXTERNAL LSAME, IZAMAX, DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, MAX, MIN
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZUNT01
      // ..
      // .. Executable Statements ..

      // Check inputs

      INFO = 0
      if ( LSAME( RC, 'R' ) ) {
         IRC = 0
      } else if ( LSAME( RC, 'C' ) ) {
         IRC = 1
      } else {
         IRC = -1
      }
      if ( IRC == -1 ) {
         INFO = -1
      } else if ( MU.LT.0 ) {
         INFO = -2
      } else if ( MV.LT.0 ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( K.LT.0 .OR. K.GT.MAX( MU, MV ) ) {
         INFO = -5
      } else if ( ( IRC == 0 .AND. LDU.LT.MAX( 1, MU ) ) .OR. ( IRC == 1 .AND. LDU.LT.MAX( 1, N ) ) ) {
         INFO = -7
      } else if ( ( IRC == 0 .AND. LDV.LT.MAX( 1, MV ) ) .OR. ( IRC == 1 .AND. LDV.LT.MAX( 1, N ) ) ) {
         INFO = -9
      }
      if ( INFO.NE.0 ) {
         xerbla('ZUNT03', -INFO );
         RETURN
      }

      // Initialize result

      RESULT = ZERO
      if (MU == 0 .OR. MV == 0 .OR. N == 0) RETURN;

      // Machine constants

      ULP = DLAMCH( 'Precision' )

      if ( IRC == 0 ) {

         // Compare rows

         RES1 = ZERO
         for (I = 1; I <= K; I++) { // 20
            LMX = IZAMAX( N, U( I, 1 ), LDU )
            if ( V( I, LMX ) == DCMPLX( ZERO ) ) {
               SV = ONE
            } else {
               SV = ABS( V( I, LMX ) ) / V( I, LMX )
            }
            if ( U( I, LMX ) == DCMPLX( ZERO ) ) {
               SU = ONE
            } else {
               SU = ABS( U( I, LMX ) ) / U( I, LMX )
            }
            S = SV / SU
            for (J = 1; J <= N; J++) { // 10
               RES1 = MAX( RES1, ABS( U( I, J )-S*V( I, J ) ) )
            } // 10
         } // 20
         RES1 = RES1 / ( DBLE( N )*ULP )

         // Compute orthogonality of rows of V.

         zunt01('Rows', MV, N, V, LDV, WORK, LWORK, RWORK, RES2 );

      } else {

         // Compare columns

         RES1 = ZERO
         for (I = 1; I <= K; I++) { // 40
            LMX = IZAMAX( N, U( 1, I ), 1 )
            if ( V( LMX, I ) == DCMPLX( ZERO ) ) {
               SV = ONE
            } else {
               SV = ABS( V( LMX, I ) ) / V( LMX, I )
            }
            if ( U( LMX, I ) == DCMPLX( ZERO ) ) {
               SU = ONE
            } else {
               SU = ABS( U( LMX, I ) ) / U( LMX, I )
            }
            S = SV / SU
            for (J = 1; J <= N; J++) { // 30
               RES1 = MAX( RES1, ABS( U( J, I )-S*V( J, I ) ) )
            } // 30
         } // 40
         RES1 = RES1 / ( DBLE( N )*ULP )

         // Compute orthogonality of columns of V.

         zunt01('Columns', N, MV, V, LDV, WORK, LWORK, RWORK, RES2 );
      }

      RESULT = MIN( MAX( RES1, RES2 ), ONE / ULP )
      RETURN

      // End of ZUNT03

      }
