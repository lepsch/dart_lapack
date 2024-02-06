      int                I1, I2;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ICHAR


      // Determine the character set.

      I1 = ICHAR( 'A' );
      I2 = ICHAR( 'a' );
      if ( I2-I1 == 32 ) {
         WRITE( *, * ) ' ASCII String    set';
      } else {
         WRITE( *, * ) ' Non-ASCII String    set, IOFF should be ',I2-I1;
      }

      // Test lsame.

      if( !lsame( 'A', 'A' ) ) WRITE( *, 9999 )'A', 'A';
      if( !lsame( 'A', 'a' ) ) WRITE( *, 9999 )'A', 'a';
      if( !lsame( 'a', 'A' ) ) WRITE( *, 9999 )'a', 'A';
      if( !lsame( 'a', 'a' ) ) WRITE( *, 9999 )'a', 'a';
      if( lsame( 'A', 'B' ) ) WRITE( *, 9998 )'A', 'B';
      if( lsame( 'A', 'b' ) ) WRITE( *, 9998 )'A', 'b';
      if( lsame( 'a', 'B' ) ) WRITE( *, 9998 )'a', 'B';
      if( lsame( 'a', 'b' ) ) WRITE( *, 9998 )'a', 'b';
      if( lsame( 'O', '/' ) ) WRITE( *, 9998 )'O', '/';
      if( lsame( '/', 'O' ) ) WRITE( *, 9998 )'/', 'O';
      if( lsame( 'o', '/' ) ) WRITE( *, 9998 )'o', '/';
      IF( lsame( '/', 'o' ) ) WRITE( *, 9998 )'/', 'o';
      WRITE( *, * )' Tests completed';

 9999 FORMAT( ' *** Error:  lsame( ', A1, , '${.a1}) is false ' );
 9998 FORMAT( ' *** Error:  lsame( ', A1, , '${.a1}) is true ' );
      }
