      RECURSIVE SUBROUTINE SGELQT3( M, N, A, LDA, T, LDT, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int       INFO, LDA, M, N, LDT;
      // ..
      // .. Array Arguments ..
      REAL      A( LDA, * ), T( LDT, * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      REAL   ONE
      PARAMETER ( ONE = 1.0E+00 )
      // ..
      // .. Local Scalars ..
      int       I, I1, J, J1, M1, M2, IINFO;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARFG, STRMM, SGEMM, XERBLA
      // ..
      // .. Executable Statements ..
*
      INFO = 0
      IF( M .LT. 0 ) THEN
         INFO = -1
      ELSE IF( N .LT. M ) THEN
         INFO = -2
      ELSE IF( LDA .LT. MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LDT .LT. MAX( 1, M ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGELQT3', -INFO )
         RETURN
      END IF
*
      IF( M.EQ.1 ) THEN
*
         // Compute Householder transform when M=1
*
         CALL SLARFG( N, A( 1, 1 ), A( 1, MIN( 2, N ) ), LDA, T( 1, 1 ) )
*
      ELSE
*
         // Otherwise, split A into blocks...
*
         M1 = M/2
         M2 = M-M1
         I1 = MIN( M1+1, M )
         J1 = MIN( M+1, N )
*
         // Compute A(1:M1,1:N) <- (Y1,R1,T1), where Q1 = I - Y1 T1 Y1^H
*
         CALL SGELQT3( M1, N, A, LDA, T, LDT, IINFO )
*
         // Compute A(J1:M,1:N) = Q1^H A(J1:M,1:N) [workspace: T(1:N1,J1:N)]
*
         DO I=1,M2
            DO J=1,M1
               T(  I+M1, J ) = A( I+M1, J )
            END DO
         END DO
         CALL STRMM( 'R', 'U', 'T', 'U', M2, M1, ONE, A, LDA, T( I1, 1 ), LDT )
*
         CALL SGEMM( 'N', 'T', M2, M1, N-M1, ONE, A( I1, I1 ), LDA, A( 1, I1 ), LDA, ONE, T( I1, 1 ), LDT)
*
         CALL STRMM( 'R', 'U', 'N', 'N', M2, M1, ONE, T, LDT, T( I1, 1 ), LDT )
*
         CALL SGEMM( 'N', 'N', M2, N-M1, M1, -ONE, T( I1, 1 ), LDT, A( 1, I1 ), LDA, ONE, A( I1, I1 ), LDA )
*
         CALL STRMM( 'R', 'U', 'N', 'U', M2, M1 , ONE, A, LDA, T( I1, 1 ), LDT )
*
         DO I=1,M2
            DO J=1,M1
               A(  I+M1, J ) = A( I+M1, J ) - T( I+M1, J )
               T( I+M1, J )=0
            END DO
         END DO
*
         // Compute A(J1:M,J1:N) <- (Y2,R2,T2) where Q2 = I - Y2 T2 Y2^H
*
         CALL SGELQT3( M2, N-M1, A( I1, I1 ), LDA, T( I1, I1 ), LDT, IINFO )
*
         // Compute T3 = T(J1:N1,1:N) = -T1 Y1^H Y2 T2
*
         DO I=1,M2
            DO J=1,M1
               T( J, I+M1  ) = (A( J, I+M1 ))
            END DO
         END DO
*
         CALL STRMM( 'R', 'U', 'T', 'U', M1, M2, ONE, A( I1, I1 ), LDA, T( 1, I1 ), LDT )
*
         CALL SGEMM( 'N', 'T', M1, M2, N-M, ONE, A( 1, J1 ), LDA, A( I1, J1 ), LDA, ONE, T( 1, I1 ), LDT )
*
         CALL STRMM( 'L', 'U', 'N', 'N', M1, M2, -ONE, T, LDT, T( 1, I1 ), LDT )
*
         CALL STRMM( 'R', 'U', 'N', 'N', M1, M2, ONE, T( I1, I1 ), LDT, T( 1, I1 ), LDT )
*
*
*
         // Y = (Y1,Y2); L = [ L1            0  ];  T = [T1 T3]
                          // [ A(1:N1,J1:N)  L2 ]       [ 0 T2]
*
      END IF
*
      RETURN
*
      // End of SGELQT3
*
      END
