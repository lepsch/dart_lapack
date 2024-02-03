      double           FUNCTION ZQRT12( M, N, A, LDA, S, WORK, LWORK, RWORK );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             RWORK( * ), S( * );
      COMPLEX*16         A( LDA, * ), WORK( LWORK )
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
      // .. Local Arrays ..
      double             DUMMY( 1 );
      // ..
      // .. External Functions ..
      double             DASUM, DLAMCH, DNRM2, ZLANGE;
      // EXTERNAL DASUM, DLAMCH, DNRM2, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DBDSQR, DLASCL, XERBLA, ZGEBD2, ZLASCL, ZLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX, MIN
      // ..
      // .. Executable Statements ..

      ZQRT12 = ZERO

      // Test that enough workspace is supplied

      if ( LWORK.LT.M*N+2*MIN( M, N )+MAX( M, N ) ) {
         xerbla('ZQRT12', 7 );
         RETURN
      }

      // Quick return if possible

      MN = MIN( M, N )
      IF( MN.LE.ZERO ) RETURN

      NRMSVL = DNRM2( MN, S, 1 )

      // Copy upper triangle of A into work

      zlaset('Full', M, N, DCMPLX( ZERO ), DCMPLX( ZERO ), WORK, M );
      DO J = 1, N
         DO I = 1, MIN( J, M )
            WORK( ( J-1 )*M+I ) = A( I, J )
         END DO
      END DO

      // Get machine parameters

      SMLNUM = DLAMCH( 'S' ) / DLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM

      // Scale work if max entry outside range [SMLNUM,BIGNUM]

      ANRM = ZLANGE( 'M', M, N, WORK, M, DUMMY )
      ISCL = 0
      if ( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         zlascl('G', 0, 0, ANRM, SMLNUM, M, N, WORK, M, INFO );
         ISCL = 1
      } else if ( ANRM.GT.BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         zlascl('G', 0, 0, ANRM, BIGNUM, M, N, WORK, M, INFO );
         ISCL = 1
      }

      if ( ANRM.NE.ZERO ) {

         // Compute SVD of work

         zgebd2(M, N, WORK, M, RWORK( 1 ), RWORK( MN+1 ), WORK( M*N+1 ), WORK( M*N+MN+1 ), WORK( M*N+2*MN+1 ), INFO )          CALL DBDSQR( 'Upper', MN, 0, 0, 0, RWORK( 1 ), RWORK( MN+1 ), DUMMY, MN, DUMMY, 1, DUMMY, MN, RWORK( 2*MN+1 ), INFO );

         if ( ISCL.EQ.1 ) {
            if ( ANRM.GT.BIGNUM ) {
               dlascl('G', 0, 0, BIGNUM, ANRM, MN, 1, RWORK( 1 ), MN, INFO );
            }
            if ( ANRM.LT.SMLNUM ) {
               dlascl('G', 0, 0, SMLNUM, ANRM, MN, 1, RWORK( 1 ), MN, INFO );
            }
         }

      } else {

         DO I = 1, MN
            RWORK( I ) = ZERO
         END DO
      }

      // Compare s and singular values of work

      daxpy(MN, -ONE, S, 1, RWORK( 1 ), 1 );
      ZQRT12 = DASUM( MN, RWORK( 1 ), 1 ) / ( DLAMCH( 'Epsilon' )*DBLE( MAX( M, N ) ) )

      IF( NRMSVL.NE.ZERO ) ZQRT12 = ZQRT12 / NRMSVL

      RETURN

      // End of ZQRT12

      }
