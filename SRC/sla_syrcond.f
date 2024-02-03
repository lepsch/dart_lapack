      REAL FUNCTION SLA_SYRCOND( UPLO, N, A, LDA, AF, LDAF, IPIV, CMODE, C, INFO, WORK, IWORK )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                N, LDA, LDAF, INFO, CMODE;
      // ..
      // .. Array Arguments
      int                IWORK( * ), IPIV( * );
      REAL               A( LDA, * ), AF( LDAF, * ), WORK( * ), C( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      String             NORMIN;
      int                KASE, I, J;
      REAL               AINVNM, SMLNUM, TMP
      bool               UP;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH
      // EXTERNAL LSAME, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACN2, XERBLA, SSYTRS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      SLA_SYRCOND = 0.0

      INFO = 0
      if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      } else if ( LDAF.LT.MAX( 1, N ) ) {
         INFO = -6
      }
      if ( INFO.NE.0 ) {
         xerbla('SLA_SYRCOND', -INFO );
         RETURN
      }
      if ( N.EQ.0 ) {
         SLA_SYRCOND = 1.0
         RETURN
      }
      UP = .FALSE.
      IF ( LSAME( UPLO, 'U' ) ) UP = .TRUE.

      // Compute the equilibration matrix R such that
      // inv(R)*A*C has unit 1-norm.

      if ( UP ) {
         for (I = 1; I <= N; I++) {
            TMP = 0.0
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
                  TMP = TMP + ABS( A( J, I ) / C( J ) )
               END DO
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + ABS( A( I, J ) / C( J ) )
               END DO
            }
            WORK( 2*N+I ) = TMP
         END DO
      } else {
         for (I = 1; I <= N; I++) {
            TMP = 0.0
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
                  TMP = TMP + ABS( A( I, J) / C( J ) )
               END DO
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + ABS( A( J, I) / C( J ) )
               END DO
            }
            WORK( 2*N+I ) = TMP
         END DO
      }

      // Estimate the norm of inv(op(A)).

      SMLNUM = SLAMCH( 'Safe minimum' )
      AINVNM = 0.0
      NORMIN = 'N'

      KASE = 0
      } // 10
      slacn2(N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE );
      if ( KASE.NE.0 ) {
         if ( KASE.EQ.2 ) {

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) * WORK( 2*N+I )
            END DO

            if ( UP ) {
               ssytrs('U', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            } else {
               ssytrs('L', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
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
               ssytrs('U', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            } else {
               ssytrs('L', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            }

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) * WORK( 2*N+I )
            END DO
         }

         GO TO 10
      }

      // Compute the estimate of the reciprocal condition number.

      IF( AINVNM .NE. 0.0 ) SLA_SYRCOND = ( 1.0 / AINVNM )

      RETURN

      // End of SLA_SYRCOND

      }
