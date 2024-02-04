      void dgehrd(N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHI, ILO, INFO, LDA, LWORK, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), TAU( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      int                NBMAX, LDT, TSIZE;
      const              NBMAX = 64, LDT = NBMAX+1, TSIZE = LDT*NBMAX ;
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IB, IINFO, IWT, J, LDWORK, LWKOPT, NB, NBMIN, NH, NX;
      double             EI;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DGEHD2, DGEMM, DLAHR2, DLARFB, DTRMM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0;
      LQUERY = ( LWORK == -1 );
      if ( N < 0 ) {
         INFO = -1;
      } else if ( ILO < 1 || ILO > max( 1, N ) ) {
         INFO = -2;
      } else if ( IHI < min( ILO, N ) || IHI > N ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      } else if ( LWORK < max( 1, N ) && !LQUERY ) {
         INFO = -8;
      }

      NH = IHI - ILO + 1;
      if ( INFO == 0 ) {

         // Compute the workspace requirements

         if ( NH <= 1 ) {
            LWKOPT = 1;
         } else {
            NB = min( NBMAX, ILAENV( 1, 'DGEHRD', ' ', N, ILO, IHI, -1 ) );
            LWKOPT = N*NB + TSIZE;
         }
         WORK[1] = LWKOPT;
      }

      if ( INFO != 0 ) {
         xerbla('DGEHRD', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Set elements 1:ILO-1 and IHI:N-1 of TAU to zero

      for (I = 1; I <= ILO - 1; I++) { // 10
         TAU[I] = ZERO;
      } // 10
      for (I = max( 1, IHI ); I <= N - 1; I++) { // 20
         TAU[I] = ZERO;
      } // 20

      // Quick return if possible

      if ( NH <= 1 ) {
         WORK[1] = 1;
         return;
      }

      // Determine the block size

      NB = min( NBMAX, ILAENV( 1, 'DGEHRD', ' ', N, ILO, IHI, -1 ) );
      NBMIN = 2;
      if ( NB > 1 && NB < NH ) {

         // Determine when to cross over from blocked to unblocked code
         // (last block is always handled by unblocked code)

         NX = max( NB, ILAENV( 3, 'DGEHRD', ' ', N, ILO, IHI, -1 ) );
         if ( NX < NH ) {

            // Determine if workspace is large enough for blocked code

            if ( LWORK < LWKOPT ) {

               // Not enough workspace to use optimal NB:  determine the
               // minimum value of NB, and reduce NB or force use of
               // unblocked code

               NBMIN = max( 2, ILAENV( 2, 'DGEHRD', ' ', N, ILO, IHI, -1 ) );
               if ( LWORK >= (N*NBMIN + TSIZE) ) {
                  NB = (LWORK-TSIZE) / N;
               } else {
                  NB = 1;
               }
            }
         }
      }
      LDWORK = N;

      if ( NB < NBMIN || NB >= NH ) {

         // Use unblocked code below

         I = ILO;

      } else {

         // Use blocked code

         IWT = 1 + N*NB;
         for (I = ILO; NB < 0 ? I >= IHI - 1 - NX : I <= IHI - 1 - NX; I += NB) { // 40
            IB = min( NB, IHI-I );

            // Reduce columns i:i+ib-1 to Hessenberg form, returning the
            // matrices V and T of the block reflector H = I - V*T*V**T
            // which performs the reduction, and also the matrix Y = A*V*T

            dlahr2(IHI, I, IB, A( 1, I ), LDA, TAU( I ), WORK( IWT ), LDT, WORK, LDWORK );

            // Apply the block reflector H to A(1:ihi,i+ib:ihi) from the
            // right, computing  A := A - Y * V**T. V(i+ib,ib-1) must be set
            // to 1

            EI = A( I+IB, I+IB-1 );
            A[I+IB, I+IB-1] = ONE;
            dgemm('No transpose', 'Transpose', IHI, IHI-I-IB+1, IB, -ONE, WORK, LDWORK, A( I+IB, I ), LDA, ONE, A( 1, I+IB ), LDA );
            A[I+IB, I+IB-1] = EI;

            // Apply the block reflector H to A(1:i,i+1:i+ib-1) from the
            // right

            dtrmm('Right', 'Lower', 'Transpose', 'Unit', I, IB-1, ONE, A( I+1, I ), LDA, WORK, LDWORK );
            for (J = 0; J <= IB-2; J++) { // 30
               daxpy(I, -ONE, WORK( LDWORK*J+1 ), 1, A( 1, I+J+1 ), 1 );
            } // 30

            // Apply the block reflector H to A(i+1:ihi,i+ib:n) from the
            // left

            dlarfb('Left', 'Transpose', 'Forward', 'Columnwise', IHI-I, N-I-IB+1, IB, A( I+1, I ), LDA, WORK( IWT ), LDT, A( I+1, I+IB ), LDA, WORK, LDWORK );
         } // 40
      }

      // Use unblocked code to reduce the rest of the matrix

      dgehd2(N, I, IHI, A, LDA, TAU, WORK, IINFO );

      WORK[1] = LWKOPT;

      return;
      }
