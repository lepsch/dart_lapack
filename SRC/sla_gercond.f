      REAL FUNCTION SLA_GERCOND( TRANS, N, A, LDA, AF, LDAF, IPIV, CMODE, C, INFO, WORK, IWORK )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                N, LDA, LDAF, INFO, CMODE;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      REAL               A( LDA, * ), AF( LDAF, * ), WORK( * ), C( * )
*    ..

*  =====================================================================

      // .. Local Scalars ..
      bool               NOTRANS;
      int                KASE, I, J;
      REAL               AINVNM, TMP
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACN2, SGETRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      SLA_GERCOND = 0.0

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
         CALL XERBLA( 'SLA_GERCOND', -INFO )
         RETURN
      }
      if ( N.EQ.0 ) {
         SLA_GERCOND = 1.0
         RETURN
      }

      // Compute the equilibration matrix R such that
      // inv(R)*A*C has unit 1-norm.

      if (NOTRANS) {
         DO I = 1, N
            TMP = 0.0
            if ( CMODE .EQ. 1 ) {
               DO J = 1, N
                  TMP = TMP + ABS( A( I, J ) * C( J ) )
               END DO
            } else if ( CMODE .EQ. 0 ) {
               DO J = 1, N
                  TMP = TMP + ABS( A( I, J ) )
               END DO
            } else {
               DO J = 1, N
                  TMP = TMP + ABS( A( I, J ) / C( J ) )
               END DO
            }
            WORK( 2*N+I ) = TMP
         END DO
      } else {
         DO I = 1, N
            TMP = 0.0
            if ( CMODE .EQ. 1 ) {
               DO J = 1, N
                  TMP = TMP + ABS( A( J, I ) * C( J ) )
               END DO
            } else if ( CMODE .EQ. 0 ) {
               DO J = 1, N
                  TMP = TMP + ABS( A( J, I ) )
               END DO
            } else {
               DO J = 1, N
                  TMP = TMP + ABS( A( J, I ) / C( J ) )
               END DO
            }
            WORK( 2*N+I ) = TMP
         END DO
      }

      // Estimate the norm of inv(op(A)).

      AINVNM = 0.0

      KASE = 0
   10 CONTINUE
      CALL SLACN2( N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE )
      if ( KASE.NE.0 ) {
         if ( KASE.EQ.2 ) {

            // Multiply by R.

            DO I = 1, N
               WORK(I) = WORK(I) * WORK(2*N+I)
            END DO

            if (NOTRANS) {
               CALL SGETRS( 'No transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO )
            } else {
               CALL SGETRS( 'Transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO )
            }

            // Multiply by inv(C).

            if ( CMODE .EQ. 1 ) {
               DO I = 1, N
                  WORK( I ) = WORK( I ) / C( I )
               END DO
            } else if ( CMODE .EQ. -1 ) {
               DO I = 1, N
                  WORK( I ) = WORK( I ) * C( I )
               END DO
            }
         } else {

            // Multiply by inv(C**T).

            if ( CMODE .EQ. 1 ) {
               DO I = 1, N
                  WORK( I ) = WORK( I ) / C( I )
               END DO
            } else if ( CMODE .EQ. -1 ) {
               DO I = 1, N
                  WORK( I ) = WORK( I ) * C( I )
               END DO
            }

            if (NOTRANS) {
               CALL SGETRS( 'Transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO )
            } else {
               CALL SGETRS( 'No transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO )
            }

            // Multiply by R.

            DO I = 1, N
               WORK( I ) = WORK( I ) * WORK( 2*N+I )
            END DO
         }
         GO TO 10
      }

      // Compute the estimate of the reciprocal condition number.

      IF( AINVNM .NE. 0.0 ) SLA_GERCOND = ( 1.0 / AINVNM )

      RETURN

      // End of SLA_GERCOND

      }
