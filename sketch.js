var setupGlobalVariables = function() {
  // set canvas size to fill the window
  xRes = windowWidth;
  yRes = windowHeight;
  minRes = min( xRes , yRes );
  maxRes = max( xRes , yRes );
  
  // size of the "field" extends past the edges of the canvas.
  xMin = 0.01*xRes;
  xMax = 0.99*xRes;
  yMin = 0.01*yRes;
  yMax = 0.99*yRes;
  
  
  // Threshold distance. Lines will be draw between particles
  // colser than this value
  distThreshold = 80;
  distMinimum = 30;
  // Distance at which lines begin to fade
  fadeThreshold = 70;
  // controls how "bright" the image is
  alphaFactor = 0.8;
  
  // R, G, B and alpha values for line color
  lineR = 255;
  lineG = 255;
  lineB = 255;
  lineAlpha = 255;
  lineColor = color( lineR , lineG , lineB , lineAlpha );
  // background color alpha
  bgAlpha = 255;
  fillAlpha = 255;
  
  // time between frames, for physics simulation
  dt = 1.0 / ( 40 );
  
  // constands for physics simulation
  edgeSpringConstant = 100;
  frictionConstant = 0.001;
  universalConstant = 250000;
  epsilon = 40;

  // number of particles - "dots"
  numDots = 150;
  
  // values for randomizing initial particle velocities
  minVel = 0*minRes;
  maxVel = 0*minRes;
  
  // mass of particles
  avgMass = 50;
  sdMass = 0;
  minMass = -50;
  maxMass = 50;

  // log of product of all distances at init
  initLogDist = 0;
  
  //
  frameCounter = 0;
  generateNew = false;
  
  clearFirstTime = true;
  startTime = 0;
  waitTime = 3000;
};

// class definition for Dots
class Dots{
  // constructor
  constructor( N ) {
    this.N = N;
    this.X = new Array( this.N );
    this.V = new Array( this.N );
    this.A = new Array( this.N );
    this.M = new Array( this.N );
    var nSquared = this.N*this.N;
    this.D = [];
    
    for( var i = 0 ; i < this.N ; i++ ) {
      this.X[i] = createVector( random(xMin,xMax) , random(yMin , yMax) );
      this.V[i] = p5.Vector.random2D();
      this.V[i].mult( random(minVel,maxVel) );
      this.A[i] = ( createVector( 0 , 0 ) );
      this.M[i] = ( random( minMass , maxMass) );
      /*
      if( random(0,1)>0.25 ) {
        this.M[i] = avgMass;
      } else {
        this.M[i] = -avgMass;
      }
      */
    }
  }
}
    
// method for Dots: updates the distance matrix
Dots.prototype.updateDistances = function() {
  for( var i = 0 ; i < this.N - 1 ; i++ ) {
    for( var j = i ; j < this.N ; j++ ) {
      var d = this.X[i].dist( this.X[j] );
      this.D[i*this.N + j] = d;
      this.D[j*this.N + 1] = d;
    }
  }
};

// method for Dots: retrieves the distance between dots i and j
Dots.prototype.getDist = function( i , j ) {
  var d = this.D[i*this.N + j];
  return d;
};

// method for Dots: returns the sum of all distances stored
Dots.prototype.totalDist = function( ) {
  var r = 1;
  for( var i = 0 ; i < this.N - 1 ; i++ ) {
    for( var j = i+1 ; j < this.N ; j++ ) {
      r += log(this.D[i*this.N + j]);
    }
  }
  return r;
};

// method for Dots: zeros out accelerations
Dots.prototype.zeroAccelerations = function() {
  for( var i = 0 ; i < this.N ; i++ ) {
    this.A[i] = createVector( 0 , 0 );
  }
};

//method for Dots: applies friction forces to all dots
Dots.prototype.applyFrictionForces = function() {
  for( var i = 0 ; i < this.N ; i++ ) {
    dA = createVector( this.V[i].x , this.V[i].y );
    dA.mult( -frictionConstant*dA.mag() );
    this.A[i].add( dA );
  }
};

// method for Dots: applies attractive/repulsive forces between all dots
Dots.prototype.applyMutualForces = function() {
  for( var i = 0 ; i < this.N - 1 ; i++ ) {
    for( var j = i + 1 ; j < this.N ; j++ ) {
      var d = this.getDist( i , j );
      var f = universalConstant / pow( d * d + epsilon * epsilon , 1.5 );
      var dAi = p5.Vector.sub( this.X[i] , this.X[j] );
      dAi.normalize();
      var dAj = createVector( dAi.x , dAi.y );
      var mj = this.M[j];
      var mi = this.M[i];
      if( mj*mi < 0 ) {
        dAi.mult( f * abs(mj) );
        dAj.mult( -f * abs(mi) );
      } else {
        dAi.mult( -f * abs(mj) );
        dAj.mult( f * abs(mi) );
      }
      
      //dAi.mult( f * this.M[j] );
      //dAj.mult( -f * this.M[i] );
      this.A[i].add( dAi );
      this.A[j].add( dAj );
    }
  }
};

// method for Dots: applies edge forces (springy bouncy)
Dots.prototype.applyEdgeForces = function() {
  for( var i = 0 ; i < this.N ; i++ ) {
    var x = this.X[i].x;
    var y = this.X[i].y;
    var m = abs(this.M[i]);
    var f;
    var dA;
    if( x < xMin ) {
        f = edgeSpringConstant * ( xMin - x );
        dA = createVector( f / m , 0 );
        this.A[i].add( dA );
    }
    if( y < yMin ) {
        f = edgeSpringConstant * ( yMin - y );
        dA = createVector( 0, f / m );
        this.A[i].add( dA );
    }
    if( x > xMax ) {
        f = edgeSpringConstant * ( x - ( xMax ) );
        dA = createVector( -f / m , 0 );
        this.A[i].add( dA );
    }
    if( y > yMax ) {
        f = edgeSpringConstant * ( y - ( yMax ) );
        dA = createVector( 0 , -f / m );
        this.A[i].add( dA );
    }
  }
};

// evolves the physics simulation one half step (half of dt)
Dots.prototype.evolveHalfStep = function() {
  this.zeroAccelerations();
  this.updateDistances();
  this.applyMutualForces();
  this.applyEdgeForces();
  this.applyFrictionForces();
  for( var i = 0 ; i < this.N ; i++ ) {
    
    this.V[i].add( p5.Vector.mult( this.A[i] , dt / 2 ) );
  }
}

// method for Dots: evolves the physics simulation n steps
Dots.prototype.evolveFullStep = function( num ) {
  for( var n = 0 ; n < num ; n++ ) {
    for( var i = 0 ; i < this.N ; i++ ) {
      this.X[i].add( p5.Vector.mult( this.V[i] , dt ) );
    }
    this.zeroAccelerations();
    this.updateDistances();
    this.applyMutualForces();
    this.applyEdgeForces();
    this.applyFrictionForces();
    for( var i = 0 ; i < this.N ; i++ ) {
      this.V[i].add( p5.Vector.mult( this.A[i] , dt ) );
    }
  }
}

// method for Dots: draws the particles (currently not used)
Dots.prototype.drawDots = function() {
  for( var i = 0 ; i < this.N ; i++ ) {
    var x = this.X[i].x;
    var y = this.X[i].y;
    var d = 2*log(abs(this.M[i]));
    if( this.M[i] > 0 ) {
      fill( 128 , 255 , 128 , fillAlpha );
    } else {
      fill( 128 , 128 , 255 , fillAlpha );
    }
    //var d = 10;
    ellipse( x , y , d , d );
    //console.log(x , y);
    
  }
}

// method for Dots: draws lines between particles
Dots.prototype.drawDistances = function(  ) {
  for( var i = 0 ; i < this.N - 1 ; i++ ) {
    for( var j = i + 1 ; j < this.N ; j++ ) {
      var d = this.getDist( i , j );
      if( d < distThreshold && d > distMinimum ) {
        var x1 = this.X[i].x;
        var y1 = this.X[i].y;
        var x2 = this.X[j].x;
        var y2 = this.X[j].y;
        if( d > distThreshold - fadeThreshold ) {
          var fadeAlpha = lineAlpha * (distThreshold - d) / fadeThreshold;
          lineColor = color( lineR , lineG , lineB , fadeAlpha*alphaFactor );
        } else {
          lineColor = color( lineR , lineG , lineB , lineAlpha*alphaFactor );
        }
        stroke( lineColor );
        line( x1 , y1 , x2 , y2 );
      }
    }
  }
}



setup = function() {
  console.log( 'Setting up global variables...' );
  setupGlobalVariables();
  createCanvas( xRes , yRes );
  console.log( 'setting up Dots object' );
  d = new Dots( numDots );
  d.evolveHalfStep();
  initLogDist = d.totalDist();
  background( 0 , 0 , 0 , 255 );
  
  startTime = millis();
  
  // show the title screen
  textAlign( CENTER );
  fill( 255 , 255 , 255 , 255 );
  textSize( 60 );
  text("VAPOR TRAIL" , 0.5*xRes , 0.5*yRes );
  textSize( 30 );
  text( "-marthematicist-" , 0.5*xRes , 0.5*yRes + 35 );
}

draw = function() {
  // if title screen wait time has not passed, do nothing
  if( millis() - startTime < waitTime ) {
    
    //return
  }
  // clear the title screen before the first frame
  if( clearFirstTime ) {
    background( 0 , 0 , 0 , 255 );
    clearFirstTime = false;
  }
  
  // sets fill color
  fill( 255 , 255 , 255 , fillAlpha );
  noStroke();
  
  if( mouseIsPressed ) {
    d.X[0] = createVector( mouseX , mouseY );
    d.M[0] = -100*maxMass;
  } else {
    d.M[0] = maxMass;
  }
  
  
  //console.log( 'draw function: evolving full step' );
  d.evolveFullStep(2);
  
  //console.log( 'drawing dots' );
  if( true ){
    background(0 , 0 , 0 , bgAlpha);
  }
  d.drawDots();
  //d.drawDistances(  );
  
  if( d.totalDist() < 0.91*initLogDist ) {
    d = new Dots( numDots );
    d.evolveHalfStep();
    initLogDist = d.totalDist();
  }

  frameCounter++;
}

// save an image of the canvas if 's' is typed
function keyTyped() {
  if( key === 's' ) {
    saveCanvas( 'canvas' , 'jpg' );
    console.log("saved");
  }
}
