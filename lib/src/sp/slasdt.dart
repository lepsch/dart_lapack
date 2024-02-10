      void slasdt(N, LVL, ND, INODE, NDIML, NDIMR, MSUB ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LVL, MSUB, N, ND;
      int                INODE( * ), NDIML( * ), NDIMR( * );
      // ..

      double               TWO;
      const              TWO = 2.0 ;
      int                I, IL, IR, LLST, MAXN, NCRNT, NLVL;
      double               TEMP;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, LOG, MAX, REAL

      // Find the number of levels on the tree.

      MAXN = max( 1, N );
      TEMP = LOG( double( MAXN ) / REAL( MSUB+1 ) ) / LOG( TWO );
      LVL = INT( TEMP ) + 1;

      I = N / 2;
      INODE[1] = I + 1;
      NDIML[1] = I;
      NDIMR[1] = N - I - 1;
      IL = 0;
      IR = 1;
      LLST = 1;
      for (NLVL = 1; NLVL <= LVL - 1; NLVL++) { // 20

         // Constructing the tree at (NLVL+1)-st level. The number of
         // nodes created on this level is LLST * 2.

         for (I = 0; I <= LLST - 1; I++) { // 10
            IL = IL + 2;
            IR = IR + 2;
            NCRNT = LLST + I;
            NDIML[IL] = NDIML( NCRNT ) / 2;
            NDIMR[IL] = NDIML( NCRNT ) - NDIML( IL ) - 1;
            INODE[IL] = INODE( NCRNT ) - NDIMR( IL ) - 1;
            NDIML[IR] = NDIMR( NCRNT ) / 2;
            NDIMR[IR] = NDIMR( NCRNT ) - NDIML( IR ) - 1;
            INODE[IR] = INODE( NCRNT ) + NDIML( IR ) + 1;
         } // 10
         LLST = LLST*2;
      } // 20
      ND = LLST*2 - 1;

      }
