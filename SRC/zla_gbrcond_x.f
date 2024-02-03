      double           FUNCTION ZLA_GBRCOND_X( TRANS, N, KL, KU, AB, LDAB, AFB, LDAFB, IPIV, X, INFO, WORK, RWORK );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                N, KL, KU, KD, KE, LDAB, LDAFB, INFO;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         AB( LDAB, * ), AFB( LDAFB, * ), WORK( * ), X( * )
      double             RWORK( * );


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
      // EXTERNAL ZLACN2, ZGBTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      ZLA_GBRCOND_X = 0.0D+0

      INFO = 0
      NOTRANS = LSAME( TRANS, 'N' )
      if ( .NOT. NOTRANS .AND. .NOT. LSAME(TRANS, 'T') .AND. .NOT. LSAME( TRANS, 'C' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( KL.LT.0 .OR. KL.GT.N-1 ) {
         INFO = -3
      } else if ( KU.LT.0 .OR. KU.GT.N-1 ) {
         INFO = -4
      } else if ( LDAB.LT.KL+KU+1 ) {
         INFO = -6
      } else if ( LDAFB.LT.2*KL+KU+1 ) {
         INFO = -8
      }
      if ( INFO.NE.0 ) {
         xerbla('ZLA_GBRCOND_X', -INFO );
         RETURN
      }

      // Compute norm of op(A)*op2(C).

      KD = KU + 1
      KE = KL + 1
      ANORM = 0.0D+0
      if ( NOTRANS ) {
         DO I = 1, N
            TMP = 0.0D+0
            DO J = MAX( I-KL, 1 ), MIN( I+KU, N )
               TMP = TMP + CABS1( AB( KD+I-J, J) * X( J ) )
            END DO
            RWORK( I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      } else {
         DO I = 1, N
            TMP = 0.0D+0
            DO J = MAX( I-KL, 1 ), MIN( I+KU, N )
               TMP = TMP + CABS1( AB( KE-I+J, I ) * X( J ) )
            END DO
            RWORK( I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      }

      // Quick return if possible.

      if ( N.EQ.0 ) {
         ZLA_GBRCOND_X = 1.0D+0
         RETURN
      } else if ( ANORM .EQ. 0.0D+0 ) {
         RETURN
      }

      // Estimate the norm of inv(op(A)).

      AINVNM = 0.0D+0

      KASE = 0
   10 CONTINUE
      zlacn2(N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE );
      if ( KASE.NE.0 ) {
         if ( KASE.EQ.2 ) {

            // Multiply by R.

            DO I = 1, N
               WORK( I ) = WORK( I ) * RWORK( I )
            END DO

            if ( NOTRANS ) {
               zgbtrs('No transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
            } else {
               zgbtrs('Conjugate transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
            ENDIF

            // Multiply by inv(X).

            DO I = 1, N
               WORK( I ) = WORK( I ) / X( I )
            END DO
         } else {

            // Multiply by inv(X**H).

            DO I = 1, N
               WORK( I ) = WORK( I ) / X( I )
            END DO

            if ( NOTRANS ) {
               zgbtrs('Conjugate transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
            } else {
               zgbtrs('No transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
            }

            // Multiply by R.

            DO I = 1, N
               WORK( I ) = WORK( I ) * RWORK( I )
            END DO
         }
         GO TO 10
      }

      // Compute the estimate of the reciprocal condition number.

      IF( AINVNM .NE. 0.0D+0 ) ZLA_GBRCOND_X = 1.0D+0 / AINVNM

      RETURN

      // End of ZLA_GBRCOND_X

      }
