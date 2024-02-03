      double           FUNCTION DLA_GERCOND( TRANS, N, A, LDA, AF, LDAF, IPIV, CMODE, C, INFO, WORK, IWORK );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                N, LDA, LDAF, INFO, CMODE;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      double             A( LDA, * ), AF( LDAF, * ), WORK( * ), C( * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               NOTRANS;
      int                KASE, I, J;
      double             AINVNM, TMP;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLACN2, DGETRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      DLA_GERCOND = 0.0D+0

      INFO = 0
      NOTRANS = LSAME( TRANS, 'N' )
      if ( .NOT. NOTRANS .AND. .NOT. LSAME(TRANS, 'T') .AND. .NOT. LSAME(TRANS, 'C') ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      } else if ( LDAF.LT.MAX( 1, N ) ) {
         INFO = -6
      }
      if ( INFO.NE.0 ) {
         xerbla('DLA_GERCOND', -INFO );
         RETURN
      }
      if ( N.EQ.0 ) {
         DLA_GERCOND = 1.0D+0
         RETURN
      }

      // Compute the equilibration matrix R such that
      // inv(R)*A*C has unit 1-norm.

      if (NOTRANS) {
         for (I = 1; I <= N; I++) {
            TMP = 0.0D+0
            if ( CMODE .EQ. 1 ) {
               for (J = 1; J <= N; J++) {
                  TMP = TMP + ABS( A( I, J ) * C( J ) )
               END DO
            } else if ( CMODE .EQ. 0 ) {
               for (J = 1; J <= N; J++) {
                  TMP = TMP + ABS( A( I, J ) )
               END DO
            } else {
               for (J = 1; J <= N; J++) {
                  TMP = TMP + ABS( A( I, J ) / C( J ) )
               END DO
            }
            WORK( 2*N+I ) = TMP
         END DO
      } else {
         for (I = 1; I <= N; I++) {
            TMP = 0.0D+0
            if ( CMODE .EQ. 1 ) {
               for (J = 1; J <= N; J++) {
                  TMP = TMP + ABS( A( J, I ) * C( J ) )
               END DO
            } else if ( CMODE .EQ. 0 ) {
               for (J = 1; J <= N; J++) {
                  TMP = TMP + ABS( A( J, I ) )
               END DO
            } else {
               for (J = 1; J <= N; J++) {
                  TMP = TMP + ABS( A( J, I ) / C( J ) )
               END DO
            }
            WORK( 2*N+I ) = TMP
         END DO
      }

      // Estimate the norm of inv(op(A)).

      AINVNM = 0.0D+0

      KASE = 0
   10 CONTINUE
      dlacn2(N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE );
      if ( KASE.NE.0 ) {
         if ( KASE.EQ.2 ) {

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK(I) = WORK(I) * WORK(2*N+I)
            END DO

            if (NOTRANS) {
               dgetrs('No transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            } else {
               dgetrs('Transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
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

            if (NOTRANS) {
               dgetrs('Transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            } else {
               dgetrs('No transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            }

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) * WORK( 2*N+I )
            END DO
         }
         GO TO 10
      }

      // Compute the estimate of the reciprocal condition number.

      IF( AINVNM .NE. 0.0D+0 ) DLA_GERCOND = ( 1.0D+0 / AINVNM )

      RETURN

      // End of DLA_GERCOND

      }
