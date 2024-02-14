      void cgelsy(final int M, final int N, final int NRHS, final Matrix<double> A_, final int LDA, final Matrix<double> B_, final int LDB, final int JPVT, final int RCOND, final int RANK, final Array<double> WORK_, final int LWORK, final Array<double> RWORK_, final Box<int> INFO,) {
  final A = A_.dim();
  final B = B_.dim();
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LDB, LWORK, M, N, NRHS, RANK;
      double               RCOND;
      int                JPVT( * );
      double               RWORK( * );
      Complex            A( LDA, * ), B( LDB, * ), WORK( * );
      // ..

      int                IMAX, IMIN;
      const              IMAX = 1, IMIN = 2 ;
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      bool               LQUERY;
      int                I, IASCL, IBSCL, ISMAX, ISMIN, J, LWKOPT, MN, NB, NB1, NB2, NB3, NB4;
      double               ANRM, BIGNUM, BNRM, SMAX, SMAXPR, SMIN, SMINPR, SMLNUM, WSIZE;
      Complex            C1, C2, S1, S2;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CGEQP3, CLAIC1, CLASCL, CLASET, CTRSM, CTZRZF, CUNMQR, CUNMRZ, XERBLA
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      //- REAL               CLANGE, SLAMCH;
      // EXTERNAL CLANGE, ILAENV, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL, CMPLX

      MN = min( M, N );
      ISMIN = MN + 1;
      ISMAX = 2*MN + 1;

      // Test the input arguments.

      INFO = 0;
      NB1 = ilaenv( 1, 'CGEQRF', ' ', M, N, -1, -1 );
      NB2 = ilaenv( 1, 'CGERQF', ' ', M, N, -1, -1 );
      NB3 = ilaenv( 1, 'CUNMQR', ' ', M, N, NRHS, -1 );
      NB4 = ilaenv( 1, 'CUNMRQ', ' ', M, N, NRHS, -1 );
      NB = max( NB1, NB2, NB3, NB4 );
      LWKOPT = max( 1, MN+2*N+NB*(N+1), 2*MN+NB*NRHS );
      WORK[1] = CMPLX( LWKOPT );
      LQUERY = ( LWORK == -1 );
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -5;
      } else if ( LDB < max( 1, M, N ) ) {
         INFO = -7;
      } else if ( LWORK < ( MN+max( 2*MN, N+1, MN+NRHS ) ) && !LQUERY ) {
         INFO = -12;
      }

      if ( INFO != 0 ) {
         xerbla('CGELSY', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( min( M, N, NRHS ) == 0 ) {
         RANK = 0;
         return;
      }

      // Get machine parameters

      SMLNUM = SLAMCH( 'S' ) / SLAMCH( 'P' );
      BIGNUM = ONE / SMLNUM;

      // Scale A, B if max entries outside range [SMLNUM,BIGNUM]

      ANRM = CLANGE( 'M', M, N, A, LDA, RWORK );
      IASCL = 0;
      if ( ANRM > ZERO && ANRM < SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         clascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO );
         IASCL = 1;
      } else if ( ANRM > BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         clascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO );
         IASCL = 2;
      } else if ( ANRM == ZERO ) {

         // Matrix all zero. Return zero solution.

         claset('F', max( M, N ), NRHS, CZERO, CZERO, B, LDB );
         RANK = 0;
         GO TO 70;
      }

      BNRM = CLANGE( 'M', M, NRHS, B, LDB, RWORK );
      IBSCL = 0;
      if ( BNRM > ZERO && BNRM < SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         clascl('G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO );
         IBSCL = 1;
      } else if ( BNRM > BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         clascl('G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO );
         IBSCL = 2;
      }

      // Compute QR factorization with column pivoting of A:
      //    A * P = Q * R

      cgeqp3(M, N, A, LDA, JPVT, WORK( 1 ), WORK( MN+1 ), LWORK-MN, RWORK, INFO );
      WSIZE = MN + double( WORK( MN+1 ) );

      // complex workspace: MN+NB*(N+1). real workspace 2*N.
      // Details of Householder rotations stored in WORK(1:MN).

      // Determine RANK using incremental condition estimation

      WORK[ISMIN] = CONE;
      WORK[ISMAX] = CONE;
      SMAX = ( A( 1, 1 ) ).abs();
      SMIN = SMAX;
      if ( ( A( 1, 1 ) ).abs() == ZERO ) {
         RANK = 0;
         claset('F', max( M, N ), NRHS, CZERO, CZERO, B, LDB );
         GO TO 70;
      } else {
         RANK = 1;
      }

      } // 10
      if ( RANK < MN ) {
         I = RANK + 1;
         claic1(IMIN, RANK, WORK( ISMIN ), SMIN, A( 1, I ), A( I, I ), SMINPR, S1, C1 );
         claic1(IMAX, RANK, WORK( ISMAX ), SMAX, A( 1, I ), A( I, I ), SMAXPR, S2, C2 );

         if ( SMAXPR*RCOND <= SMINPR ) {
            for (I = 1; I <= RANK; I++) { // 20
               WORK[ISMIN+I-1] = S1*WORK( ISMIN+I-1 );
               WORK[ISMAX+I-1] = S2*WORK( ISMAX+I-1 );
            } // 20
            WORK[ISMIN+RANK] = C1;
            WORK[ISMAX+RANK] = C2;
            SMIN = SMINPR;
            SMAX = SMAXPR;
            RANK = RANK + 1;
            GO TO 10;
         }
      }

      // complex workspace: 3*MN.

      // Logically partition R = [ R11 R12 ]
      //                         [  0  R22 ]
      // where R11 = R(1:RANK,1:RANK)

      // [R11,R12] = [ T11, 0 ] * Y

      if (RANK < N) ctzrzf( RANK, N, A, LDA, WORK( MN+1 ), WORK( 2*MN+1 ), LWORK-2*MN, INFO );

      // complex workspace: 2*MN.
      // Details of Householder rotations stored in WORK(MN+1:2*MN)

      // B(1:M,1:NRHS) := Q**H * B(1:M,1:NRHS)

      cunmqr('Left', 'Conjugate transpose', M, NRHS, MN, A, LDA, WORK( 1 ), B, LDB, WORK( 2*MN+1 ), LWORK-2*MN, INFO );
      WSIZE = max( WSIZE, 2*MN+double( WORK( 2*MN+1 ) ) );

      // complex workspace: 2*MN+NB*NRHS.

      // B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS)

      ctrsm('Left', 'Upper', 'No transpose', 'Non-unit', RANK, NRHS, CONE, A, LDA, B, LDB );

      for (J = 1; J <= NRHS; J++) { // 40
         for (I = RANK + 1; I <= N; I++) { // 30
            B[I][J] = CZERO;
         } // 30
      } // 40

      // B(1:N,1:NRHS) := Y**H * B(1:N,1:NRHS)

      if ( RANK < N ) {
         cunmrz('Left', 'Conjugate transpose', N, NRHS, RANK, N-RANK, A, LDA, WORK( MN+1 ), B, LDB, WORK( 2*MN+1 ), LWORK-2*MN, INFO );
      }

      // complex workspace: 2*MN+NRHS.

      // B(1:N,1:NRHS) := P * B(1:N,1:NRHS)

      for (J = 1; J <= NRHS; J++) { // 60
         for (I = 1; I <= N; I++) { // 50
            WORK[JPVT( I )] = B( I, J );
         } // 50
         ccopy(N, WORK( 1 ), 1, B( 1, J ), 1 );
      } // 60

      // complex workspace: N.

      // Undo scaling

      if ( IASCL == 1 ) {
         clascl('G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO );
         clascl('U', 0, 0, SMLNUM, ANRM, RANK, RANK, A, LDA, INFO );
      } else if ( IASCL == 2 ) {
         clascl('G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO );
         clascl('U', 0, 0, BIGNUM, ANRM, RANK, RANK, A, LDA, INFO );
      }
      if ( IBSCL == 1 ) {
         clascl('G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO );
      } else if ( IBSCL == 2 ) {
         clascl('G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO );
      }

      } // 70
      WORK[1] = CMPLX( LWKOPT );

      }
