      SUBROUTINE CGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                INFO, LDA, LDB, LWORK, M, N, NRHS;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CZERO
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ) )
      // ..
      // .. Local Scalars ..
      bool               LQUERY, TPSD;
      int                BROW, I, IASCL, IBSCL, J, MN, NB, SCLLEN, WSIZE;
      REAL               ANRM, BIGNUM, BNRM, SMLNUM
      // ..
      // .. Local Arrays ..
      REAL               RWORK( 1 )
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               CLANGE, SLAMCH
      // EXTERNAL LSAME, ILAENV, CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGELQF, CGEQRF, CLASCL, CLASET, CTRTRS, CUNMLQ, CUNMQR, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      // Test the input arguments.

      INFO = 0
      MN = MIN( M, N )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.( LSAME( TRANS, 'N' ) .OR. LSAME( TRANS, 'C' ) ) ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, M, N ) ) THEN
         INFO = -8
      ELSE IF( LWORK.LT.MAX( 1, MN+MAX( MN, NRHS ) ) .AND. .NOT.LQUERY ) THEN
         INFO = -10
      END IF

      // Figure out optimal block size

      IF( INFO.EQ.0 .OR. INFO.EQ.-10 ) THEN

         TPSD = .TRUE.
         IF( LSAME( TRANS, 'N' ) ) TPSD = .FALSE.

         IF( M.GE.N ) THEN
            NB = ILAENV( 1, 'CGEQRF', ' ', M, N, -1, -1 )
            IF( TPSD ) THEN
               NB = MAX( NB, ILAENV( 1, 'CUNMQR', 'LN', M, NRHS, N, -1 ) )
            ELSE
               NB = MAX( NB, ILAENV( 1, 'CUNMQR', 'LC', M, NRHS, N, -1 ) )
            END IF
         ELSE
            NB = ILAENV( 1, 'CGELQF', ' ', M, N, -1, -1 )
            IF( TPSD ) THEN
               NB = MAX( NB, ILAENV( 1, 'CUNMLQ', 'LC', N, NRHS, M, -1 ) )
            ELSE
               NB = MAX( NB, ILAENV( 1, 'CUNMLQ', 'LN', N, NRHS, M, -1 ) )
            END IF
         END IF

         WSIZE = MAX( 1, MN + MAX( MN, NRHS )*NB )
         WORK( 1 ) = REAL( WSIZE )

      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGELS ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Quick return if possible

      IF( MIN( M, N, NRHS ).EQ.0 ) THEN
         CALL CLASET( 'Full', MAX( M, N ), NRHS, CZERO, CZERO, B, LDB )
         RETURN
      END IF

      // Get machine parameters

      SMLNUM = SLAMCH( 'S' ) / SLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM

      // Scale A, B if max element outside range [SMLNUM,BIGNUM]

      ANRM = CLANGE( 'M', M, N, A, LDA, RWORK )
      IASCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN

         // Scale matrix norm up to SMLNUM

         CALL CLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
         IASCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN

         // Scale matrix norm down to BIGNUM

         CALL CLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
         IASCL = 2
      ELSE IF( ANRM.EQ.ZERO ) THEN

         // Matrix all zero. Return zero solution.

         CALL CLASET( 'F', MAX( M, N ), NRHS, CZERO, CZERO, B, LDB )
         GO TO 50
      END IF

      BROW = M
      IF( TPSD ) BROW = N
      BNRM = CLANGE( 'M', BROW, NRHS, B, LDB, RWORK )
      IBSCL = 0
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN

         // Scale matrix norm up to SMLNUM

         CALL CLASCL( 'G', 0, 0, BNRM, SMLNUM, BROW, NRHS, B, LDB, INFO )
         IBSCL = 1
      ELSE IF( BNRM.GT.BIGNUM ) THEN

         // Scale matrix norm down to BIGNUM

         CALL CLASCL( 'G', 0, 0, BNRM, BIGNUM, BROW, NRHS, B, LDB, INFO )
         IBSCL = 2
      END IF

      IF( M.GE.N ) THEN

         // compute QR factorization of A

         CALL CGEQRF( M, N, A, LDA, WORK( 1 ), WORK( MN+1 ), LWORK-MN, INFO )

         // workspace at least N, optimally N*NB

         IF( .NOT.TPSD ) THEN

            // Least-Squares Problem min || A * X - B ||

            // B(1:M,1:NRHS) := Q**H * B(1:M,1:NRHS)

            CALL CUNMQR( 'Left', 'Conjugate transpose', M, NRHS, N, A, LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO )

            // workspace at least NRHS, optimally NRHS*NB

            // B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)

            CALL CTRTRS( 'Upper', 'No transpose', 'Non-unit', N, NRHS, A, LDA, B, LDB, INFO )

            IF( INFO.GT.0 ) THEN
               RETURN
            END IF

            SCLLEN = N

         ELSE

            // Underdetermined system of equations A**T * X = B

            // B(1:N,1:NRHS) := inv(R**H) * B(1:N,1:NRHS)

            CALL CTRTRS( 'Upper', 'Conjugate transpose','Non-unit', N, NRHS, A, LDA, B, LDB, INFO )

            IF( INFO.GT.0 ) THEN
               RETURN
            END IF

            // B(N+1:M,1:NRHS) = ZERO

            DO 20 J = 1, NRHS
               DO 10 I = N + 1, M
                  B( I, J ) = CZERO
   10          CONTINUE
   20       CONTINUE

            // B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS)

            CALL CUNMQR( 'Left', 'No transpose', M, NRHS, N, A, LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO )

            // workspace at least NRHS, optimally NRHS*NB

            SCLLEN = M

         END IF

      ELSE

         // Compute LQ factorization of A

         CALL CGELQF( M, N, A, LDA, WORK( 1 ), WORK( MN+1 ), LWORK-MN, INFO )

         // workspace at least M, optimally M*NB.

         IF( .NOT.TPSD ) THEN

            // underdetermined system of equations A * X = B

            // B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)

            CALL CTRTRS( 'Lower', 'No transpose', 'Non-unit', M, NRHS, A, LDA, B, LDB, INFO )

            IF( INFO.GT.0 ) THEN
               RETURN
            END IF

            // B(M+1:N,1:NRHS) = 0

            DO 40 J = 1, NRHS
               DO 30 I = M + 1, N
                  B( I, J ) = CZERO
   30          CONTINUE
   40       CONTINUE

            // B(1:N,1:NRHS) := Q(1:N,:)**H * B(1:M,1:NRHS)

            CALL CUNMLQ( 'Left', 'Conjugate transpose', N, NRHS, M, A, LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO )

            // workspace at least NRHS, optimally NRHS*NB

            SCLLEN = N

         ELSE

            // overdetermined system min || A**H * X - B ||

            // B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)

            CALL CUNMLQ( 'Left', 'No transpose', N, NRHS, M, A, LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO )

            // workspace at least NRHS, optimally NRHS*NB

            // B(1:M,1:NRHS) := inv(L**H) * B(1:M,1:NRHS)

            CALL CTRTRS( 'Lower', 'Conjugate transpose', 'Non-unit', M, NRHS, A, LDA, B, LDB, INFO )

            IF( INFO.GT.0 ) THEN
               RETURN
            END IF

            SCLLEN = M

         END IF

      END IF

      // Undo scaling

      IF( IASCL.EQ.1 ) THEN
         CALL CLASCL( 'G', 0, 0, ANRM, SMLNUM, SCLLEN, NRHS, B, LDB, INFO )
      ELSE IF( IASCL.EQ.2 ) THEN
         CALL CLASCL( 'G', 0, 0, ANRM, BIGNUM, SCLLEN, NRHS, B, LDB, INFO )
      END IF
      IF( IBSCL.EQ.1 ) THEN
         CALL CLASCL( 'G', 0, 0, SMLNUM, BNRM, SCLLEN, NRHS, B, LDB, INFO )
      ELSE IF( IBSCL.EQ.2 ) THEN
         CALL CLASCL( 'G', 0, 0, BIGNUM, BNRM, SCLLEN, NRHS, B, LDB, INFO )
      END IF

   50 CONTINUE
      WORK( 1 ) = REAL( WSIZE )

      RETURN

      // End of CGELS

      END
