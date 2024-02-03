      double           FUNCTION ZLA_PORCOND_C( UPLO, N, A, LDA, AF, LDAF, C, CAPPLY, INFO, WORK, RWORK );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      bool               CAPPLY;
      int                N, LDA, LDAF, INFO;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), AF( LDAF, * ), WORK( * )
      double             C( * ), RWORK( * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                KASE;
      double             AINVNM, ANORM, TMP;
      int                I, J;
      bool               UP, UPPER;
      COMPLEX*16         ZDUM
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLACN2, ZPOTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, REAL, DIMAG
      // ..
      // .. Statement Functions ..
      double           CABS1;
      // ..
      // .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      ZLA_PORCOND_C = 0.0D+0

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      } else if ( LDAF.LT.MAX( 1, N ) ) {
         INFO = -6
      }
      if ( INFO.NE.0 ) {
         xerbla('ZLA_PORCOND_C', -INFO );
         RETURN
      }
      UP = .FALSE.
      IF ( LSAME( UPLO, 'U' ) ) UP = .TRUE.

      // Compute norm of op(A)*op2(C).

      ANORM = 0.0D+0
      if ( UP ) {
         for (I = 1; I <= N; I++) {
            TMP = 0.0D+0
            if ( CAPPLY ) {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + CABS1( A( J, I ) ) / C( J )
               END DO
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + CABS1( A( I, J ) ) / C( J )
               END DO
            } else {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + CABS1( A( J, I ) )
               END DO
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + CABS1( A( I, J ) )
               END DO
            }
            RWORK( I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      } else {
         for (I = 1; I <= N; I++) {
            TMP = 0.0D+0
            if ( CAPPLY ) {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + CABS1( A( I, J ) ) / C( J )
               END DO
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + CABS1( A( J, I ) ) / C( J )
               END DO
            } else {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + CABS1( A( I, J ) )
               END DO
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + CABS1( A( J, I ) )
               END DO
            }
            RWORK( I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      }

      // Quick return if possible.

      if ( N.EQ.0 ) {
         ZLA_PORCOND_C = 1.0D+0
         RETURN
      } else if ( ANORM .EQ. 0.0D+0 ) {
         RETURN
      }

      // Estimate the norm of inv(op(A)).

      AINVNM = 0.0D+0

      KASE = 0
      } // 10
      zlacn2(N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE );
      if ( KASE.NE.0 ) {
         if ( KASE.EQ.2 ) {

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) * RWORK( I )
            END DO

            if ( UP ) {
               zpotrs('U', N, 1, AF, LDAF, WORK, N, INFO );
            } else {
               zpotrs('L', N, 1, AF, LDAF, WORK, N, INFO );
            }

            // Multiply by inv(C).

            if ( CAPPLY ) {
               for (I = 1; I <= N; I++) {
                  WORK( I ) = WORK( I ) * C( I )
               END DO
            }
         } else {

            // Multiply by inv(C**H).

            if ( CAPPLY ) {
               for (I = 1; I <= N; I++) {
                  WORK( I ) = WORK( I ) * C( I )
               END DO
            }

            if ( UP ) {
               zpotrs('U', N, 1, AF, LDAF, WORK, N, INFO );
            } else {
               zpotrs('L', N, 1, AF, LDAF, WORK, N, INFO );
            }

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) * RWORK( I )
            END DO
         }
         GO TO 10
      }

      // Compute the estimate of the reciprocal condition number.

      IF( AINVNM .NE. 0.0D+0 ) ZLA_PORCOND_C = 1.0D+0 / AINVNM

      RETURN

      // End of ZLA_PORCOND_C

      }
