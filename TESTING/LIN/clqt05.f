      SUBROUTINE CLQT05(M,N,L,NB,RESULT)
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int     LWORK, M, N, L, NB, LDT;
      // .. Return values ..
      REAL RESULT(6)
*
*  =====================================================================
*
      // ..
      // .. Local allocatable arrays
      COMPLEX, ALLOCATABLE :: AF(:,:), Q(:,:), R(:,:), WORK( : ), T(:,:), CF(:,:), DF(:,:), A(:,:), C(:,:), D(:,:)
      REAL, ALLOCATABLE :: RWORK(:)
*
      // .. Parameters ..
      REAL ZERO
      COMPLEX       ONE, CZERO
      PARAMETER( ZERO = 0.0, ONE = (1.0,0.0), CZERO=(0.0,0.0) )
      // ..
      // .. Local Scalars ..
      int     INFO, J, K, N2, NP1,i;
      REAL   ANORM, EPS, RESID, CNORM, DNORM
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 );
      // ..
      // .. External Functions ..
      REAL SLAMCH
      REAL CLANGE, CLANSY
      bool     LSAME;
      // EXTERNAL SLAMCH, CLANGE, CLANSY, LSAME
      // ..
      // .. Data statements ..
      DATA ISEED / 1988, 1989, 1990, 1991 /
*
      EPS = SLAMCH( 'Epsilon' )
      K = M
      N2 = M+N
      IF( N.GT.0 ) THEN
         NP1 = M+1
      ELSE
         NP1 = 1
      END IF
      LWORK = N2*N2*NB
*
      // Dynamically allocate all arrays
*
      ALLOCATE(A(M,N2),AF(M,N2),Q(N2,N2),R(N2,N2),RWORK(N2), WORK(LWORK),T(NB,M),C(N2,M),CF(N2,M), D(M,N2),DF(M,N2) )
*
      // Put random stuff into A
*
      LDT=NB
      CALL CLASET( 'Full', M, N2, CZERO, CZERO, A, M )
      CALL CLASET( 'Full', NB, M, CZERO, CZERO, T, NB )
      DO J=1,M
         CALL CLARNV( 2, ISEED, M-J+1, A( J, J ) )
      END DO
      IF( N.GT.0 ) THEN
         DO J=1,N-L
            CALL CLARNV( 2, ISEED, M, A( 1, MIN(N+M,M+1) + J - 1 ) )
         END DO
      END IF
      IF( L.GT.0 ) THEN
         DO J=1,L
            CALL CLARNV( 2, ISEED, M-J+1, A( J, MIN(N+M,N+M-L+1) + J - 1 ) )
         END DO
      END IF
*
      // Copy the matrix A to the array AF.
*
      CALL CLACPY( 'Full', M, N2, A, M, AF, M )
*
      // Factor the matrix A in the array AF.
*
      CALL CTPLQT( M,N,L,NB,AF,M,AF(1,NP1),M,T,LDT,WORK,INFO)
*
      // Generate the (M+N)-by-(M+N) matrix Q by applying H to I
*
      CALL CLASET( 'Full', N2, N2, CZERO, ONE, Q, N2 )
      CALL CGEMLQT( 'L', 'N', N2, N2, K, NB, AF, M, T, LDT, Q, N2, WORK, INFO )
*
      // Copy L
*
      CALL CLASET( 'Full', N2, N2, CZERO, CZERO, R, N2 )
      CALL CLACPY( 'Lower', M, N2, AF, M, R, N2 )
*
      // Compute |L - A*Q*C| / |A| and store in RESULT(1)
*
      CALL CGEMM( 'N', 'C', M, N2, N2, -ONE,  A, M, Q, N2, ONE, R, N2)
      ANORM = CLANGE( '1', M, N2, A, M, RWORK )
      RESID = CLANGE( '1', M, N2, R, N2, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = RESID / (EPS*ANORM*MAX(1,N2))
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
      // Compute |I - Q*Q'| and store in RESULT(2)
*
      CALL CLASET( 'Full', N2, N2, CZERO, ONE, R, N2 )
      CALL CHERK( 'U', 'N', N2, N2, REAL(-ONE), Q, N2, REAL(ONE), R, N2 )
      RESID = CLANSY( '1', 'Upper', N2, R, N2, RWORK )
      RESULT( 2 ) = RESID / (EPS*MAX(1,N2))
*
      // Generate random m-by-n matrix C and a copy CF
*
      CALL CLASET( 'Full', N2, M, CZERO, ONE, C, N2 )
      DO J=1,M
         CALL CLARNV( 2, ISEED, N2, C( 1, J ) )
      END DO
      CNORM = CLANGE( '1', N2, M, C, N2, RWORK)
      CALL CLACPY( 'Full', N2, M, C, N2, CF, N2 )
*
      // Apply Q to C as Q*C
*
      CALL CTPMLQT( 'L','N', N,M,K,L,NB,AF(1, NP1),M,T,LDT,CF,N2, CF(NP1,1),N2,WORK,INFO)
*
      // Compute |Q*C - Q*C| / |C|
*
      CALL CGEMM( 'N', 'N', N2, M, N2, -ONE, Q, N2, C, N2, ONE, CF, N2 )
      RESID = CLANGE( '1', N2, M, CF, N2, RWORK )
      IF( CNORM.GT.ZERO ) THEN
         RESULT( 3 ) = RESID / (EPS*MAX(1,N2)*CNORM)
      ELSE
         RESULT( 3 ) = ZERO
      END IF

*
      // Copy C into CF again
*
      CALL CLACPY( 'Full', N2, M, C, N2, CF, N2 )
*
      // Apply Q to C as QT*C
*
      CALL CTPMLQT( 'L','C',N,M,K,L,NB,AF(1,NP1),M,T,LDT,CF,N2, CF(NP1,1),N2,WORK,INFO)
*
      // Compute |QT*C - QT*C| / |C|
*
      CALL CGEMM('C','N',N2,M,N2,-ONE,Q,N2,C,N2,ONE,CF,N2)
      RESID = CLANGE( '1', N2, M, CF, N2, RWORK )

      IF( CNORM.GT.ZERO ) THEN
         RESULT( 4 ) = RESID / (EPS*MAX(1,N2)*CNORM)
      ELSE
         RESULT( 4 ) = ZERO
      END IF
*
      // Generate random m-by-n matrix D and a copy DF
*
      DO J=1,N2
         CALL CLARNV( 2, ISEED, M, D( 1, J ) )
      END DO
      DNORM = CLANGE( '1', M, N2, D, M, RWORK)
      CALL CLACPY( 'Full', M, N2, D, M, DF, M )
*
      // Apply Q to D as D*Q
*
      CALL CTPMLQT('R','N',M,N,K,L,NB,AF(1,NP1),M,T,LDT,DF,M, DF(1,NP1),M,WORK,INFO)
*
      // Compute |D*Q - D*Q| / |D|
*
      CALL CGEMM('N','N',M,N2,N2,-ONE,D,M,Q,N2,ONE,DF,M)
      RESID = CLANGE('1',M, N2,DF,M,RWORK )
      IF( CNORM.GT.ZERO ) THEN
         RESULT( 5 ) = RESID / (EPS*MAX(1,N2)*DNORM)
      ELSE
         RESULT( 5 ) = ZERO
      END IF
*
      // Copy D into DF again
*
      CALL CLACPY('Full',M,N2,D,M,DF,M )
*
      // Apply Q to D as D*QT
*
      CALL CTPMLQT('R','C',M,N,K,L,NB,AF(1,NP1),M,T,LDT,DF,M, DF(1,NP1),M,WORK,INFO)

*
      // Compute |D*QT - D*QT| / |D|
*
      CALL CGEMM( 'N', 'C', M, N2, N2, -ONE, D, M, Q, N2, ONE, DF, M )
      RESID = CLANGE( '1', M, N2, DF, M, RWORK )
      IF( CNORM.GT.ZERO ) THEN
         RESULT( 6 ) = RESID / (EPS*MAX(1,N2)*DNORM)
      ELSE
         RESULT( 6 ) = ZERO
      END IF
*
      // Deallocate all arrays
*
      DEALLOCATE ( A, AF, Q, R, RWORK, WORK, T, C, D, CF, DF)
      RETURN
      END
