      double           FUNCTION DLA_PORCOND( UPLO, N, A, LDA, AF, LDAF, CMODE, C, INFO, WORK, IWORK );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                N, LDA, LDAF, INFO, CMODE;
      double             A( LDA, * ), AF( LDAF, * ), WORK( * ), C( * );
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                KASE, I, J;
      double             AINVNM, TMP;
      bool               UP;
      // ..
      // .. Array Arguments ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLACN2, DPOTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      DLA_PORCOND = 0.0D+0

      INFO = 0
      if ( N.LT.0 ) {
         INFO = -2
      }
      if ( INFO.NE.0 ) {
         xerbla('DLA_PORCOND', -INFO );
         RETURN
      }

      if ( N.EQ.0 ) {
         DLA_PORCOND = 1.0D+0
         RETURN
      }
      UP = .FALSE.
      IF ( LSAME( UPLO, 'U' ) ) UP = .TRUE.

      // Compute the equilibration matrix R such that
      // inv(R)*A*C has unit 1-norm.

      if ( UP ) {
         for (I = 1; I <= N; I++) {
            TMP = 0.0D+0
            if ( CMODE .EQ. 1 ) {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + ABS( A( J, I ) * C( J ) )
               END DO
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + ABS( A( I, J ) * C( J ) )
               END DO
            } else if ( CMODE .EQ. 0 ) {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + ABS( A( J, I ) )
               END DO
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + ABS( A( I, J ) )
               END DO
            } else {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + ABS( A( J ,I ) / C( J ) )
               END DO
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + ABS( A( I, J ) / C( J ) )
               END DO
            }
            WORK( 2*N+I ) = TMP
         END DO
      } else {
         for (I = 1; I <= N; I++) {
            TMP = 0.0D+0
            if ( CMODE .EQ. 1 ) {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + ABS( A( I, J ) * C( J ) )
               END DO
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + ABS( A( J, I ) * C( J ) )
               END DO
            } else if ( CMODE .EQ. 0 ) {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + ABS( A( I, J ) )
               END DO
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + ABS( A( J, I ) )
               END DO
            } else {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + ABS( A( I, J ) / C( J ) )
               END DO
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + ABS( A( J, I ) / C( J ) )
               END DO
            }
            WORK( 2*N+I ) = TMP
         END DO
      }

      // Estimate the norm of inv(op(A)).

      AINVNM = 0.0D+0

      KASE = 0
      } // 10
      dlacn2(N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE );
      if ( KASE.NE.0 ) {
         if ( KASE.EQ.2 ) {

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) * WORK( 2*N+I )
            END DO

            if (UP) {
               dpotrs('Upper', N, 1, AF, LDAF, WORK, N, INFO );
            } else {
               dpotrs('Lower', N, 1, AF, LDAF, WORK, N, INFO );
            }

            // Multiply by inv(C).

            if ( CMODE .EQ. 1 ) {
               for (I = 1; I <= N; I++) {
                  WORK( I ) = WORK( I ) / C( I )
               END DO
            } else if ( CMODE .EQ. -1 ) {
               for (I = 1; I <= N; I++) {
                  WORK( I ) = WORK( I ) * C( I )
               END DO
            }
         } else {

            // Multiply by inv(C**T).

            if ( CMODE .EQ. 1 ) {
               for (I = 1; I <= N; I++) {
                  WORK( I ) = WORK( I ) / C( I )
               END DO
            } else if ( CMODE .EQ. -1 ) {
               for (I = 1; I <= N; I++) {
                  WORK( I ) = WORK( I ) * C( I )
               END DO
            }

            if ( UP ) {
               dpotrs('Upper', N, 1, AF, LDAF, WORK, N, INFO );
            } else {
               dpotrs('Lower', N, 1, AF, LDAF, WORK, N, INFO );
            }

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) * WORK( 2*N+I )
            END DO
         }
         GO TO 10
      }

      // Compute the estimate of the reciprocal condition number.

      IF( AINVNM .NE. 0.0D+0 ) DLA_PORCOND = ( 1.0D+0 / AINVNM )

      RETURN

      // End of DLA_PORCOND

      }
