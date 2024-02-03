      double           FUNCTION ZLA_GERCOND_C( TRANS, N, A, LDA, AF, LDAF, IPIV, C, CAPPLY, INFO, WORK, RWORK );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      bool               CAPPLY;
      int                N, LDA, LDAF, INFO;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         A( LDA, * ), AF( LDAF, * ), WORK( * )
      double             C( * ), RWORK( * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               NOTRANS;
      int                KASE, I, J;
      double             AINVNM, ANORM, TMP;
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
      // EXTERNAL ZLACN2, ZGETRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, REAL, DIMAG
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..
      ZLA_GERCOND_C = 0.0D+0

      INFO = 0
      NOTRANS = LSAME( TRANS, 'N' )
      if ( .NOT. NOTRANS .AND. .NOT. LSAME( TRANS, 'T' ) .AND. .NOT. LSAME( TRANS, 'C' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      } else if ( LDAF.LT.MAX( 1, N ) ) {
         INFO = -6
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'ZLA_GERCOND_C', -INFO )
         RETURN
      }

      // Compute norm of op(A)*op2(C).

      ANORM = 0.0D+0
      if ( NOTRANS ) {
         DO I = 1, N
            TMP = 0.0D+0
            if ( CAPPLY ) {
               DO J = 1, N
                  TMP = TMP + CABS1( A( I, J ) ) / C( J )
               END DO
            } else {
               DO J = 1, N
                  TMP = TMP + CABS1( A( I, J ) )
               END DO
            }
            RWORK( I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      } else {
         DO I = 1, N
            TMP = 0.0D+0
            if ( CAPPLY ) {
               DO J = 1, N
                  TMP = TMP + CABS1( A( J, I ) ) / C( J )
               END DO
            } else {
               DO J = 1, N
                  TMP = TMP + CABS1( A( J, I ) )
               END DO
            }
            RWORK( I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      }

      // Quick return if possible.

      if ( N.EQ.0 ) {
         ZLA_GERCOND_C = 1.0D+0
         RETURN
      } else if ( ANORM .EQ. 0.0D+0 ) {
         RETURN
      }

      // Estimate the norm of inv(op(A)).

      AINVNM = 0.0D+0

      KASE = 0
   10 CONTINUE
      CALL ZLACN2( N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE )
      if ( KASE.NE.0 ) {
         if ( KASE.EQ.2 ) {

            // Multiply by R.

            DO I = 1, N
               WORK( I ) = WORK( I ) * RWORK( I )
            END DO

            if (NOTRANS) {
               CALL ZGETRS( 'No transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO )
            } else {
               CALL ZGETRS( 'Conjugate transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO )
            ENDIF

            // Multiply by inv(C).

            if ( CAPPLY ) {
               DO I = 1, N
                  WORK( I ) = WORK( I ) * C( I )
               END DO
            }
         } else {

            // Multiply by inv(C**H).

            if ( CAPPLY ) {
               DO I = 1, N
                  WORK( I ) = WORK( I ) * C( I )
               END DO
            }

            if ( NOTRANS ) {
               CALL ZGETRS( 'Conjugate transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO )
            } else {
               CALL ZGETRS( 'No transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO )
            }

            // Multiply by R.

            DO I = 1, N
               WORK( I ) = WORK( I ) * RWORK( I )
            END DO
         }
         GO TO 10
      }

      // Compute the estimate of the reciprocal condition number.

      IF( AINVNM .NE. 0.0D+0 ) ZLA_GERCOND_C = 1.0D+0 / AINVNM

      RETURN

      // End of ZLA_GERCOND_C

      }
