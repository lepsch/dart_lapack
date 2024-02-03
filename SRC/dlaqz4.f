      SUBROUTINE DLAQZ4( ILSCHUR, ILQ, ILZ, N, ILO, IHI, NSHIFTS, NBLOCK_DESIRED, SR, SI, SS, A, LDA, B, LDB, Q, LDQ, Z, LDZ, QC, LDQC, ZC, LDZC, WORK, LWORK, INFO )
      IMPLICIT NONE

      // Function arguments
      bool   , INTENT( IN ) :: ILSCHUR, ILQ, ILZ;
      int    , INTENT( IN ) :: N, ILO, IHI, LDA, LDB, LDQ, LDZ, LWORK, NSHIFTS, NBLOCK_DESIRED, LDQC, LDZC;
       double          , INTENT( INOUT ) :: A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * ), QC( LDQC, * ), ZC( LDZC, * ), WORK( * ), SR( * ), SI( * ), SS( * );

      int    , INTENT( OUT ) :: INFO;

      // Parameters
      double           :: ZERO, ONE, HALF;
      const    ZERO = 0.0D0, ONE = 1.0D0, HALF = 0.5D0 ;

      // Local scalars
      int     :: I, J, NS, ISTARTM, ISTOPM, SHEIGHT, SWIDTH, K, NP, ISTARTB, ISTOPB, ISHIFT, NBLOCK, NPOS;
      double           :: TEMP, V( 3 ), C1, S1, C2, S2, SWAP;

      // External functions
      // EXTERNAL :: XERBLA, DGEMM, DLAQZ1, DLAQZ2, DLASET, DLARTG, DROT, DLACPY

      INFO = 0
      if ( NBLOCK_DESIRED < NSHIFTS+1 ) {
         INFO = -8
      }
      if ( LWORK == -1 ) {
         // workspace query, quick return
         WORK( 1 ) = N*NBLOCK_DESIRED
         RETURN
      } else if ( LWORK < N*NBLOCK_DESIRED ) {
         INFO = -25
      }

      if ( INFO != 0 ) {
         xerbla('DLAQZ4', -INFO );
         RETURN
      }

      // Executable statements

      if ( NSHIFTS < 2 ) {
         RETURN
      }

      if ( ILO >= IHI ) {
         RETURN
      }

      if ( ILSCHUR ) {
         ISTARTM = 1
         ISTOPM = N
      } else {
         ISTARTM = ILO
         ISTOPM = IHI
      }

      // Shuffle shifts into pairs of real shifts and pairs
      // of complex conjugate shifts assuming complex
      // conjugate shifts are already adjacent to one
      // another

      DO I = 1, NSHIFTS-2, 2
         if ( SI( I ) != -SI( I+1 ) ) {

            SWAP = SR( I )
            SR( I ) = SR( I+1 )
            SR( I+1 ) = SR( I+2 )
            SR( I+2 ) = SWAP

            SWAP = SI( I )
            SI( I ) = SI( I+1 )
            SI( I+1 ) = SI( I+2 )
            SI( I+2 ) = SWAP

            SWAP = SS( I )
            SS( I ) = SS( I+1 )
            SS( I+1 ) = SS( I+2 )
            SS( I+2 ) = SWAP
         }
      }

      // NSHFTS is supposed to be even, but if it is odd,
      // then simply reduce it by one.  The shuffle above
      // ensures that the dropped shift is real and that
      // the remaining shifts are paired.

      NS = NSHIFTS-MOD( NSHIFTS, 2 )
      NPOS = MAX( NBLOCK_DESIRED-NS, 1 )

      // The following block introduces the shifts and chases
      // them down one by one just enough to make space for
      // the other shifts. The near-the-diagonal block is
      // of size (ns+1) x ns.

      dlaset('FULL', NS+1, NS+1, ZERO, ONE, QC, LDQC );
      dlaset('FULL', NS, NS, ZERO, ONE, ZC, LDZC );

      DO I = 1, NS, 2
         // Introduce the shift
         dlaqz1(A( ILO, ILO ), LDA, B( ILO, ILO ), LDB, SR( I ), SR( I+1 ), SI( I ), SS( I ), SS( I+1 ), V );

         TEMP = V( 2 )
         dlartg(TEMP, V( 3 ), C1, S1, V( 2 ) );
         dlartg(V( 1 ), V( 2 ), C2, S2, TEMP );
          drot(NS, A( ILO+1, ILO ), LDA, A( ILO+2, ILO ), LDA, C1, S1 );
         drot(NS, A( ILO, ILO ), LDA, A( ILO+1, ILO ), LDA, C2, S2 );
         drot(NS, B( ILO+1, ILO ), LDB, B( ILO+2, ILO ), LDB, C1, S1 );
         drot(NS, B( ILO, ILO ), LDB, B( ILO+1, ILO ), LDB, C2, S2 );
         drot(NS+1, QC( 1, 2 ), 1, QC( 1, 3 ), 1, C1, S1 );
         drot(NS+1, QC( 1, 1 ), 1, QC( 1, 2 ), 1, C2, S2 );

         // Chase the shift down
         for (J = 1; J <= NS-1-I; J++) {
             dlaqz2( true , true , J, 1, NS, IHI-ILO+1, A( ILO, ILO ), LDA, B( ILO, ILO ), LDB, NS+1, 1, QC, LDQC, NS, 1, ZC, LDZC );

         }

      }

      // Update the rest of the pencil

      // Update A(ilo:ilo+ns,ilo+ns:istopm) and B(ilo:ilo+ns,ilo+ns:istopm)
      // from the left with Qc(1:ns+1,1:ns+1)'
      SHEIGHT = NS+1
      SWIDTH = ISTOPM-( ILO+NS )+1
      if ( SWIDTH > 0 ) {
         dgemm('T', 'N', SHEIGHT, SWIDTH, SHEIGHT, ONE, QC, LDQC, A( ILO, ILO+NS ), LDA, ZERO, WORK, SHEIGHT );
         dlacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( ILO, ILO+NS ), LDA );
         dgemm('T', 'N', SHEIGHT, SWIDTH, SHEIGHT, ONE, QC, LDQC, B( ILO, ILO+NS ), LDB, ZERO, WORK, SHEIGHT );
         dlacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( ILO, ILO+NS ), LDB );
      }
      if ( ILQ ) {
         dgemm('N', 'N', N, SHEIGHT, SHEIGHT, ONE, Q( 1, ILO ), LDQ, QC, LDQC, ZERO, WORK, N );
         dlacpy('ALL', N, SHEIGHT, WORK, N, Q( 1, ILO ), LDQ );
      }

      // Update A(istartm:ilo-1,ilo:ilo+ns-1) and B(istartm:ilo-1,ilo:ilo+ns-1)
      // from the right with Zc(1:ns,1:ns)
      SHEIGHT = ILO-1-ISTARTM+1
      SWIDTH = NS
      if ( SHEIGHT > 0 ) {
         dgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, ONE, A( ISTARTM, ILO ), LDA, ZC, LDZC, ZERO, WORK, SHEIGHT );
         dlacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( ISTARTM, ILO ), LDA );
         dgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, ONE, B( ISTARTM, ILO ), LDB, ZC, LDZC, ZERO, WORK, SHEIGHT );
         dlacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( ISTARTM, ILO ), LDB );
      }
      if ( ILZ ) {
         dgemm('N', 'N', N, SWIDTH, SWIDTH, ONE, Z( 1, ILO ), LDZ, ZC, LDZC, ZERO, WORK, N );
         dlacpy('ALL', N, SWIDTH, WORK, N, Z( 1, ILO ), LDZ );
      }

      // The following block chases the shifts down to the bottom
      // right block. If possible, a shift is moved down npos
      // positions at a time

      K = ILO
      DO WHILE ( K < IHI-NS )
         NP = MIN( IHI-NS-K, NPOS )
         // Size of the near-the-diagonal block
         NBLOCK = NS+NP
         // istartb points to the first row we will be updating
         ISTARTB = K+1
         // istopb points to the last column we will be updating
         ISTOPB = K+NBLOCK-1

         dlaset('FULL', NS+NP, NS+NP, ZERO, ONE, QC, LDQC );
         dlaset('FULL', NS+NP, NS+NP, ZERO, ONE, ZC, LDZC );

         // Near the diagonal shift chase
         DO I = NS-1, 0, -2
            for (J = 0; J <= NP-1; J++) {
               // Move down the block with index k+i+j-1, updating
               // the (ns+np x ns+np) block:
               // (k:k+ns+np,k:k+ns+np-1)
               dlaqz2( true , true , K+I+J-1, ISTARTB, ISTOPB, IHI, A, LDA, B, LDB, NBLOCK, K+1, QC, LDQC, NBLOCK, K, ZC, LDZC );
            }
         }

         // Update rest of the pencil

         // Update A(k+1:k+ns+np, k+ns+np:istopm) and
         // B(k+1:k+ns+np, k+ns+np:istopm)
         // from the left with Qc(1:ns+np,1:ns+np)'
         SHEIGHT = NS+NP
         SWIDTH = ISTOPM-( K+NS+NP )+1
         if ( SWIDTH > 0 ) {
            dgemm('T', 'N', SHEIGHT, SWIDTH, SHEIGHT, ONE, QC, LDQC, A( K+1, K+NS+NP ), LDA, ZERO, WORK, SHEIGHT );
            dlacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( K+1, K+NS+NP ), LDA );
            dgemm('T', 'N', SHEIGHT, SWIDTH, SHEIGHT, ONE, QC, LDQC, B( K+1, K+NS+NP ), LDB, ZERO, WORK, SHEIGHT );
            dlacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( K+1, K+NS+NP ), LDB );
         }
         if ( ILQ ) {
            dgemm('N', 'N', N, NBLOCK, NBLOCK, ONE, Q( 1, K+1 ), LDQ, QC, LDQC, ZERO, WORK, N );
            dlacpy('ALL', N, NBLOCK, WORK, N, Q( 1, K+1 ), LDQ );
         }

         // Update A(istartm:k,k:k+ns+npos-1) and B(istartm:k,k:k+ns+npos-1)
         // from the right with Zc(1:ns+np,1:ns+np)
         SHEIGHT = K-ISTARTM+1
         SWIDTH = NBLOCK
         if ( SHEIGHT > 0 ) {
            dgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, ONE, A( ISTARTM, K ), LDA, ZC, LDZC, ZERO, WORK, SHEIGHT );
            dlacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( ISTARTM, K ), LDA );
            dgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, ONE, B( ISTARTM, K ), LDB, ZC, LDZC, ZERO, WORK, SHEIGHT );
            dlacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( ISTARTM, K ), LDB );
         }
         if ( ILZ ) {
            dgemm('N', 'N', N, NBLOCK, NBLOCK, ONE, Z( 1, K ), LDZ, ZC, LDZC, ZERO, WORK, N );
            dlacpy('ALL', N, NBLOCK, WORK, N, Z( 1, K ), LDZ );
         }

         K = K+NP

      }

      // The following block removes the shifts from the bottom right corner
      // one by one. Updates are initially applied to A(ihi-ns+1:ihi,ihi-ns:ihi).

      dlaset('FULL', NS, NS, ZERO, ONE, QC, LDQC );
      dlaset('FULL', NS+1, NS+1, ZERO, ONE, ZC, LDZC );

      // istartb points to the first row we will be updating
      ISTARTB = IHI-NS+1
      // istopb points to the last column we will be updating
      ISTOPB = IHI

      DO I = 1, NS, 2
         // Chase the shift down to the bottom right corner
         for (ISHIFT = IHI-I-1; ISHIFT <= IHI-2; ISHIFT++) {
            dlaqz2( true , true , ISHIFT, ISTARTB, ISTOPB, IHI, A, LDA, B, LDB, NS, IHI-NS+1, QC, LDQC, NS+1, IHI-NS, ZC, LDZC );
         }

      }

      // Update rest of the pencil

      // Update A(ihi-ns+1:ihi, ihi+1:istopm)
      // from the left with Qc(1:ns,1:ns)'
      SHEIGHT = NS
      SWIDTH = ISTOPM-( IHI+1 )+1
      if ( SWIDTH > 0 ) {
         dgemm('T', 'N', SHEIGHT, SWIDTH, SHEIGHT, ONE, QC, LDQC, A( IHI-NS+1, IHI+1 ), LDA, ZERO, WORK, SHEIGHT );
         dlacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( IHI-NS+1, IHI+1 ), LDA );
         dgemm('T', 'N', SHEIGHT, SWIDTH, SHEIGHT, ONE, QC, LDQC, B( IHI-NS+1, IHI+1 ), LDB, ZERO, WORK, SHEIGHT );
         dlacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( IHI-NS+1, IHI+1 ), LDB );
      }
      if ( ILQ ) {
         dgemm('N', 'N', N, NS, NS, ONE, Q( 1, IHI-NS+1 ), LDQ, QC, LDQC, ZERO, WORK, N );
         dlacpy('ALL', N, NS, WORK, N, Q( 1, IHI-NS+1 ), LDQ );
      }

      // Update A(istartm:ihi-ns,ihi-ns:ihi)
      // from the right with Zc(1:ns+1,1:ns+1)
      SHEIGHT = IHI-NS-ISTARTM+1
      SWIDTH = NS+1
      if ( SHEIGHT > 0 ) {
         dgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, ONE, A( ISTARTM, IHI-NS ), LDA, ZC, LDZC, ZERO, WORK, SHEIGHT );
         dlacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( ISTARTM, IHI-NS ), LDA );
         dgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, ONE, B( ISTARTM, IHI-NS ), LDB, ZC, LDZC, ZERO, WORK, SHEIGHT );
         dlacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( ISTARTM, IHI-NS ), LDB );
      }
      if ( ILZ ) {
         dgemm('N', 'N', N, NS+1, NS+1, ONE, Z( 1, IHI-NS ), LDZ, ZC, LDZC, ZERO, WORK, N );
         dlacpy('ALL', N, NS+1, WORK, N, Z( 1, IHI-NS ), LDZ );
      }

      END SUBROUTINE
