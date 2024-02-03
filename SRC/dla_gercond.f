      double           FUNCTION DLA_GERCOND( TRANS, N, A, LDA, AF, LDAF, IPIV, CMODE, C, INFO, WORK, IWORK );
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             TRANS;
      int                N, LDA, LDAF, INFO, CMODE;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      double             A( LDA, * ), AF( LDAF, * ), WORK( * ), C( * );
      // ..
*
*  =====================================================================
*
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
*
      DLA_GERCOND = 0.0D+0
*
      INFO = 0
      NOTRANS = LSAME( TRANS, 'N' )
      IF ( .NOT. NOTRANS .AND. .NOT. LSAME(TRANS, 'T') .AND. .NOT. LSAME(TRANS, 'C') ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDAF.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLA_GERCOND', -INFO )
         RETURN
      END IF
      IF( N.EQ.0 ) THEN
         DLA_GERCOND = 1.0D+0
         RETURN
      END IF
*
      // Compute the equilibration matrix R such that
      // inv(R)*A*C has unit 1-norm.
*
      IF (NOTRANS) THEN
         DO I = 1, N
            TMP = 0.0D+0
            IF ( CMODE .EQ. 1 ) THEN
               DO J = 1, N
                  TMP = TMP + ABS( A( I, J ) * C( J ) )
               END DO
            ELSE IF ( CMODE .EQ. 0 ) THEN
               DO J = 1, N
                  TMP = TMP + ABS( A( I, J ) )
               END DO
            ELSE
               DO J = 1, N
                  TMP = TMP + ABS( A( I, J ) / C( J ) )
               END DO
            END IF
            WORK( 2*N+I ) = TMP
         END DO
      ELSE
         DO I = 1, N
            TMP = 0.0D+0
            IF ( CMODE .EQ. 1 ) THEN
               DO J = 1, N
                  TMP = TMP + ABS( A( J, I ) * C( J ) )
               END DO
            ELSE IF ( CMODE .EQ. 0 ) THEN
               DO J = 1, N
                  TMP = TMP + ABS( A( J, I ) )
               END DO
            ELSE
               DO J = 1, N
                  TMP = TMP + ABS( A( J, I ) / C( J ) )
               END DO
            END IF
            WORK( 2*N+I ) = TMP
         END DO
      END IF
*
      // Estimate the norm of inv(op(A)).
*
      AINVNM = 0.0D+0

      KASE = 0
   10 CONTINUE
      CALL DLACN2( N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
         IF( KASE.EQ.2 ) THEN
*
            // Multiply by R.
*
            DO I = 1, N
               WORK(I) = WORK(I) * WORK(2*N+I)
            END DO

            IF (NOTRANS) THEN
               CALL DGETRS( 'No transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO )
            ELSE
               CALL DGETRS( 'Transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO )
            END IF
*
            // Multiply by inv(C).
*
            IF ( CMODE .EQ. 1 ) THEN
               DO I = 1, N
                  WORK( I ) = WORK( I ) / C( I )
               END DO
            ELSE IF ( CMODE .EQ. -1 ) THEN
               DO I = 1, N
                  WORK( I ) = WORK( I ) * C( I )
               END DO
            END IF
         ELSE
*
            // Multiply by inv(C**T).
*
            IF ( CMODE .EQ. 1 ) THEN
               DO I = 1, N
                  WORK( I ) = WORK( I ) / C( I )
               END DO
            ELSE IF ( CMODE .EQ. -1 ) THEN
               DO I = 1, N
                  WORK( I ) = WORK( I ) * C( I )
               END DO
            END IF

            IF (NOTRANS) THEN
               CALL DGETRS( 'Transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO )
            ELSE
               CALL DGETRS( 'No transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO )
            END IF
*
            // Multiply by R.
*
            DO I = 1, N
               WORK( I ) = WORK( I ) * WORK( 2*N+I )
            END DO
         END IF
         GO TO 10
      END IF
*
      // Compute the estimate of the reciprocal condition number.
*
      IF( AINVNM .NE. 0.0D+0 ) DLA_GERCOND = ( 1.0D+0 / AINVNM )
*
      RETURN
*
      // End of DLA_GERCOND
*
      END
