      double           FUNCTION DQRT12( M, N, A, LDA, S, WORK, LWORK );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), S( * ), WORK( LWORK );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, ISCL, J, MN;
      double             ANRM, BIGNUM, NRMSVL, SMLNUM;
      // ..
      // .. External Functions ..
      double             DASUM, DLAMCH, DLANGE, DNRM2;
      // EXTERNAL DASUM, DLAMCH, DLANGE, DNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DBDSQR, DGEBD2, DLASCL, DLASET, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Local Arrays ..
      double             DUMMY( 1 );
      // ..
      // .. Executable Statements ..

      DQRT12 = ZERO

      // Test that enough workspace is supplied

      if ( LWORK < MAX( M*N+4*MIN( M, N )+MAX( M, N ), M*N+2*MIN( M, N )+4*N) ) {
         xerbla('DQRT12', 7 );
         RETURN
      }

      // Quick return if possible

      MN = MIN( M, N )
      if (MN.LE.ZERO) RETURN;

      NRMSVL = DNRM2( MN, S, 1 )

      // Copy upper triangle of A into work

      dlaset('Full', M, N, ZERO, ZERO, WORK, M );
      for (J = 1; J <= N; J++) {
         DO I = 1, MIN( J, M )
            WORK( ( J-1 )*M+I ) = A( I, J )
         }
      }

      // Get machine parameters

      SMLNUM = DLAMCH( 'S' ) / DLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM

      // Scale work if max entry outside range [SMLNUM,BIGNUM]

      ANRM = DLANGE( 'M', M, N, WORK, M, DUMMY )
      ISCL = 0
      if ( ANRM > ZERO && ANRM < SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         dlascl('G', 0, 0, ANRM, SMLNUM, M, N, WORK, M, INFO );
         ISCL = 1
      } else if ( ANRM > BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         dlascl('G', 0, 0, ANRM, BIGNUM, M, N, WORK, M, INFO );
         ISCL = 1
      }

      if ( ANRM != ZERO ) {

         // Compute SVD of work

         dgebd2(M, N, WORK, M, WORK( M*N+1 ), WORK( M*N+MN+1 ), WORK( M*N+2*MN+1 ), WORK( M*N+3*MN+1 ), WORK( M*N+4*MN+1 ), INFO );
         dbdsqr('Upper', MN, 0, 0, 0, WORK( M*N+1 ), WORK( M*N+MN+1 ), DUMMY, MN, DUMMY, 1, DUMMY, MN, WORK( M*N+2*MN+1 ), INFO );

         if ( ISCL == 1 ) {
            if ( ANRM > BIGNUM ) {
               dlascl('G', 0, 0, BIGNUM, ANRM, MN, 1, WORK( M*N+1 ), MN, INFO );
            }
            if ( ANRM < SMLNUM ) {
               dlascl('G', 0, 0, SMLNUM, ANRM, MN, 1, WORK( M*N+1 ), MN, INFO );
            }
         }

      } else {

         for (I = 1; I <= MN; I++) {
            WORK( M*N+I ) = ZERO
         }
      }

      // Compare s and singular values of work

      daxpy(MN, -ONE, S, 1, WORK( M*N+1 ), 1 );

      DQRT12 = DASUM( MN, WORK( M*N+1 ), 1 ) / ( DLAMCH('Epsilon') * DBLE( MAX( M, N ) ) )

      if (NRMSVL != ZERO) DQRT12 = DQRT12 / NRMSVL;

      RETURN

      // End of DQRT12

      }
