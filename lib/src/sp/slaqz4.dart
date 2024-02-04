      void slaqz4(ILSCHUR, ILQ, ILZ, N, ILO, IHI, NSHIFTS, NBLOCK_DESIRED, SR, SI, SS, A, LDA, B, LDB, Q, LDQ, Z, LDZ, QC, LDQC, ZC, LDZC, WORK, LWORK, INFO ) {
      // IMPLICIT NONE

      // Function arguments
      bool   , INTENT( IN ) :: ILSCHUR, ILQ, ILZ;
      int    , INTENT( IN ) :: N, ILO, IHI, LDA, LDB, LDQ, LDZ, LWORK, NSHIFTS, NBLOCK_DESIRED, LDQC, LDZC;
       double, INTENT( INOUT ) :: A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * ), QC( LDQC, * ), ZC( LDZC, * ), WORK( * ), SR( * ), SI( * ), SS( * );

      int    , INTENT( OUT ) :: INFO;

      // Parameters
      double :: ZERO, ONE, HALF;
      const    ZERO = 0.0, ONE = 1.0, HALF = 0.5 ;

      // Local scalars
      int     :: I, J, NS, ISTARTM, ISTOPM, SHEIGHT, SWIDTH, K, NP, ISTARTB, ISTOPB, ISHIFT, NBLOCK, NPOS;
      double :: TEMP, V( 3 ), C1, S1, C2, S2, SWAP;

      // External functions
      // EXTERNAL :: XERBLA, SGEMM, SLAQZ1, SLAQZ2, SLASET, SLARTG, SROT, SLACPY
      double, EXTERNAL :: SROUNDUP_LWORK;

      INFO = 0;
      if ( NBLOCK_DESIRED < NSHIFTS+1 ) {
         INFO = -8;
      }
      if ( LWORK == -1 ) {
         // workspace query, quick return;
         WORK[1] = SROUNDUP_LWORK(N*NBLOCK_DESIRED);
         return;
      } else if ( LWORK < N*NBLOCK_DESIRED ) {
         INFO = -25;
      }

      if ( INFO != 0 ) {
         xerbla('SLAQZ4', -INFO );
         return;
      }

      // Executable statements

      if ( NSHIFTS < 2 ) {
         return;
      }

      if ( ILO >= IHI ) {
         return;
      }

      if ( ILSCHUR ) {
         ISTARTM = 1;
         ISTOPM = N;
      } else {
         ISTARTM = ILO;
         ISTOPM = IHI;
      }

      // Shuffle shifts into pairs of real shifts and pairs
      // of complex conjugate shifts assuming complex
      // conjugate shifts are already adjacent to one
      // another

      for (I = 1; 2 < 0 ? I >= NSHIFTS-2 : I <= NSHIFTS-2; I += 2) {
         if ( SI( I ) != -SI( I+1 ) ) {

            SWAP = SR( I );
            SR[I] = SR( I+1 );
            SR[I+1] = SR( I+2 );
            SR[I+2] = SWAP;

            SWAP = SI( I );
            SI[I] = SI( I+1 );
            SI[I+1] = SI( I+2 );
            SI[I+2] = SWAP;

            SWAP = SS( I );
            SS[I] = SS( I+1 );
            SS[I+1] = SS( I+2 );
            SS[I+2] = SWAP;
         }
      }

      // NSHFTS is supposed to be even, but if it is odd,
      // then simply reduce it by one.  The shuffle above
      // ensures that the dropped shift is real and that
      // the remaining shifts are paired.

      NS = NSHIFTS-(NSHIFTS % 2);
      NPOS = max( NBLOCK_DESIRED-NS, 1 );

      // The following block introduces the shifts and chases
      // them down one by one just enough to make space for
      // the other shifts. The near-the-diagonal block is
      // of size (ns+1) x ns.

      slaset('FULL', NS+1, NS+1, ZERO, ONE, QC, LDQC );
      slaset('FULL', NS, NS, ZERO, ONE, ZC, LDZC );

      for (I = 1; I <= NS; I += 2) {
         // Introduce the shift
         slaqz1(A( ILO, ILO ), LDA, B( ILO, ILO ), LDB, SR( I ), SR( I+1 ), SI( I ), SS( I ), SS( I+1 ), V );

         TEMP = V( 2 );
         slartg(TEMP, V( 3 ), C1, S1, V( 2 ) );
         slartg(V( 1 ), V( 2 ), C2, S2, TEMP );
          srot(NS, A( ILO+1, ILO ), LDA, A( ILO+2, ILO ), LDA, C1, S1 );
         srot(NS, A( ILO, ILO ), LDA, A( ILO+1, ILO ), LDA, C2, S2 );
         srot(NS, B( ILO+1, ILO ), LDB, B( ILO+2, ILO ), LDB, C1, S1 );
         srot(NS, B( ILO, ILO ), LDB, B( ILO+1, ILO ), LDB, C2, S2 );
         srot(NS+1, QC( 1, 2 ), 1, QC( 1, 3 ), 1, C1, S1 );
         srot(NS+1, QC( 1, 1 ), 1, QC( 1, 2 ), 1, C2, S2 );

         // Chase the shift down
         for (J = 1; J <= NS-1-I; J++) {
             slaqz2( true , true , J, 1, NS, IHI-ILO+1, A( ILO, ILO ), LDA, B( ILO, ILO ), LDB, NS+1, 1, QC, LDQC, NS, 1, ZC, LDZC );

         }

      }

      // Update the rest of the pencil

      // Update A(ilo:ilo+ns,ilo+ns:istopm) and B(ilo:ilo+ns,ilo+ns:istopm)
      // from the left with Qc(1:ns+1,1:ns+1)'
      SHEIGHT = NS+1;
      SWIDTH = ISTOPM-( ILO+NS )+1;
      if ( SWIDTH > 0 ) {
         sgemm('T', 'N', SHEIGHT, SWIDTH, SHEIGHT, ONE, QC, LDQC, A( ILO, ILO+NS ), LDA, ZERO, WORK, SHEIGHT );
         slacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( ILO, ILO+NS ), LDA );
         sgemm('T', 'N', SHEIGHT, SWIDTH, SHEIGHT, ONE, QC, LDQC, B( ILO, ILO+NS ), LDB, ZERO, WORK, SHEIGHT );
         slacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( ILO, ILO+NS ), LDB );
      }
      if ( ILQ ) {
         sgemm('N', 'N', N, SHEIGHT, SHEIGHT, ONE, Q( 1, ILO ), LDQ, QC, LDQC, ZERO, WORK, N );
         slacpy('ALL', N, SHEIGHT, WORK, N, Q( 1, ILO ), LDQ );
      }

      // Update A(istartm:ilo-1,ilo:ilo+ns-1) and B(istartm:ilo-1,ilo:ilo+ns-1)
      // from the right with Zc(1:ns,1:ns)
      SHEIGHT = ILO-1-ISTARTM+1;
      SWIDTH = NS;
      if ( SHEIGHT > 0 ) {
         sgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, ONE, A( ISTARTM, ILO ), LDA, ZC, LDZC, ZERO, WORK, SHEIGHT );
         slacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( ISTARTM, ILO ), LDA );
         sgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, ONE, B( ISTARTM, ILO ), LDB, ZC, LDZC, ZERO, WORK, SHEIGHT );
         slacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( ISTARTM, ILO ), LDB );
      }
      if ( ILZ ) {
         sgemm('N', 'N', N, SWIDTH, SWIDTH, ONE, Z( 1, ILO ), LDZ, ZC, LDZC, ZERO, WORK, N );
         slacpy('ALL', N, SWIDTH, WORK, N, Z( 1, ILO ), LDZ );
      }

      // The following block chases the shifts down to the bottom
      // right block. If possible, a shift is moved down npos
      // positions at a time

      K = ILO;
      while (K < IHI-NS) {
         NP = min( IHI-NS-K, NPOS );
         // Size of the near-the-diagonal block
         NBLOCK = NS+NP;
         // istartb points to the first row we will be updating
         ISTARTB = K+1;
         // istopb points to the last column we will be updating
         ISTOPB = K+NBLOCK-1;

         slaset('FULL', NS+NP, NS+NP, ZERO, ONE, QC, LDQC );
         slaset('FULL', NS+NP, NS+NP, ZERO, ONE, ZC, LDZC );

         // Near the diagonal shift chase
         for (I = NS-1; I >= 0; I -= 2) {
            for (J = 0; J <= NP-1; J++) {
               // Move down the block with index k+i+j-1, updating
               // the (ns+np x ns+np) block:
               // (k:k+ns+np,k:k+ns+np-1)
               slaqz2( true , true , K+I+J-1, ISTARTB, ISTOPB, IHI, A, LDA, B, LDB, NBLOCK, K+1, QC, LDQC, NBLOCK, K, ZC, LDZC );
            }
         }

         // Update rest of the pencil

         // Update A(k+1:k+ns+np, k+ns+np:istopm) and
         // B(k+1:k+ns+np, k+ns+np:istopm)
         // from the left with Qc(1:ns+np,1:ns+np)'
         SHEIGHT = NS+NP;
         SWIDTH = ISTOPM-( K+NS+NP )+1;
         if ( SWIDTH > 0 ) {
            sgemm('T', 'N', SHEIGHT, SWIDTH, SHEIGHT, ONE, QC, LDQC, A( K+1, K+NS+NP ), LDA, ZERO, WORK, SHEIGHT );
            slacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( K+1, K+NS+NP ), LDA );
            sgemm('T', 'N', SHEIGHT, SWIDTH, SHEIGHT, ONE, QC, LDQC, B( K+1, K+NS+NP ), LDB, ZERO, WORK, SHEIGHT );
            slacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( K+1, K+NS+NP ), LDB );
         }
         if ( ILQ ) {
            sgemm('N', 'N', N, NBLOCK, NBLOCK, ONE, Q( 1, K+1 ), LDQ, QC, LDQC, ZERO, WORK, N );
            slacpy('ALL', N, NBLOCK, WORK, N, Q( 1, K+1 ), LDQ );
         }

         // Update A(istartm:k,k:k+ns+npos-1) and B(istartm:k,k:k+ns+npos-1)
         // from the right with Zc(1:ns+np,1:ns+np)
         SHEIGHT = K-ISTARTM+1;
         SWIDTH = NBLOCK;
         if ( SHEIGHT > 0 ) {
            sgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, ONE, A( ISTARTM, K ), LDA, ZC, LDZC, ZERO, WORK, SHEIGHT );
            slacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( ISTARTM, K ), LDA );
            sgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, ONE, B( ISTARTM, K ), LDB, ZC, LDZC, ZERO, WORK, SHEIGHT );
            slacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( ISTARTM, K ), LDB );
         }
         if ( ILZ ) {
            sgemm('N', 'N', N, NBLOCK, NBLOCK, ONE, Z( 1, K ), LDZ, ZC, LDZC, ZERO, WORK, N );
            slacpy('ALL', N, NBLOCK, WORK, N, Z( 1, K ), LDZ );
         }

         K = K+NP;

      }

      // The following block removes the shifts from the bottom right corner
      // one by one. Updates are initially applied to A(ihi-ns+1:ihi,ihi-ns:ihi).

      slaset('FULL', NS, NS, ZERO, ONE, QC, LDQC );
      slaset('FULL', NS+1, NS+1, ZERO, ONE, ZC, LDZC );

      // istartb points to the first row we will be updating
      ISTARTB = IHI-NS+1;
      // istopb points to the last column we will be updating
      ISTOPB = IHI;

      for (I = 1; I <= NS; I += 2) {
         // Chase the shift down to the bottom right corner
         for (ISHIFT = IHI-I-1; ISHIFT <= IHI-2; ISHIFT++) {
            slaqz2( true , true , ISHIFT, ISTARTB, ISTOPB, IHI, A, LDA, B, LDB, NS, IHI-NS+1, QC, LDQC, NS+1, IHI-NS, ZC, LDZC );
         }

      }

      // Update rest of the pencil

      // Update A(ihi-ns+1:ihi, ihi+1:istopm)
      // from the left with Qc(1:ns,1:ns)'
      SHEIGHT = NS;
      SWIDTH = ISTOPM-( IHI+1 )+1;
      if ( SWIDTH > 0 ) {
         sgemm('T', 'N', SHEIGHT, SWIDTH, SHEIGHT, ONE, QC, LDQC, A( IHI-NS+1, IHI+1 ), LDA, ZERO, WORK, SHEIGHT );
         slacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( IHI-NS+1, IHI+1 ), LDA );
         sgemm('T', 'N', SHEIGHT, SWIDTH, SHEIGHT, ONE, QC, LDQC, B( IHI-NS+1, IHI+1 ), LDB, ZERO, WORK, SHEIGHT );
         slacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( IHI-NS+1, IHI+1 ), LDB );
      }
      if ( ILQ ) {
         sgemm('N', 'N', N, NS, NS, ONE, Q( 1, IHI-NS+1 ), LDQ, QC, LDQC, ZERO, WORK, N );
         slacpy('ALL', N, NS, WORK, N, Q( 1, IHI-NS+1 ), LDQ );
      }

      // Update A(istartm:ihi-ns,ihi-ns:ihi)
      // from the right with Zc(1:ns+1,1:ns+1)
      SHEIGHT = IHI-NS-ISTARTM+1;
      SWIDTH = NS+1;
      if ( SHEIGHT > 0 ) {
         sgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, ONE, A( ISTARTM, IHI-NS ), LDA, ZC, LDZC, ZERO, WORK, SHEIGHT );
         slacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( ISTARTM, IHI-NS ), LDA );
         sgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, ONE, B( ISTARTM, IHI-NS ), LDB, ZC, LDZC, ZERO, WORK, SHEIGHT );
         slacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( ISTARTM, IHI-NS ), LDB );
      }
      if ( ILZ ) {
      sgemm('N', 'N', N, NS+1, NS+1, ONE, Z( 1, IHI-NS ), LDZ, ZC, LDZC, ZERO, WORK, N );
         slacpy('ALL', N, NS+1, WORK, N, Z( 1, IHI-NS ), LDZ );
      }

      END SUBROUTINE;