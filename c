/*
  Particle filters implementation for a simple robot.

  Your goal here is to implement a simple particle filter
  algorithm for robot localization.

  A map file in .ppm format is read from disk, the map
  contains empty spaces and walls.

  A simple robot is randomly placed on this map, and
  your task is to write a particle filter algorithm to
  allow the robot to determine its location with
  high certainty.

  You must complete all sections marked

  TO DO:

  NOTE: 2 keyboard controls are provided:

  q -> quit (exit program)
  r -> reset particle set during simulation

  Written by F.J.E. for CSC C85, May 2012. Updated Aug. 15, 2014

  This robot model inspired by Sebastian Thrun's
  model in CS373.
*/

#include "ParticleFilters.h"

/**********************************************************
 GLOBAL DATA
**********************************************************/
unsigned char *map;		// Input map
unsigned char *map_b;		// Temporary frame
struct particle *robot;		// Robot
struct particle *list;		// Particle list
int sx,sy;			// Size of the map image
char name[1024];		// Name of the map
int n_particles;		// Number of particles
int windowID;               	// Glut window ID (for display)
int Win[2];                 	// window (x,y) size
int RESETflag;			// RESET particles

/**********************************************************
 PROGRAM CODE
**********************************************************/
int main(int argc, char *argv[])
{
 /*
   Main function. Usage for this program:

   ParticleFilters map_name n_particles

   Where:
    map_name is the name of a .ppm file containing the map. The map
             should be BLACK on empty (free) space, and coloured
             wherever there are obstacles or walls. Anythin not
             black is an obstacle.

    n_particles is the number of particles to simulate in [100, 50000]

   Main loads the map image, initializes a robot at a random location
    in the map, and sets up the OpenGL stuff before entering the
    filtering loop.
 */

 if (argc!=3)
 {
  fprintf(stderr,"Wrong number of parameters. Usage: ParticleFilters map_name n_particles.\n");
  exit(0);
 }

 strcpy(&name[0],argv[1]);
 n_particles=atoi(argv[2]);

 if (n_particles<100||n_particles>50000)
 {
  fprintf(stderr,"Number of particles must be in [100, 50000]\n");
  exit(0);
 }

 fprintf(stderr,"Reading input map\n");
 map=readPPMimage(name,&sx, &sy);
 if (map==NULL)
 {
  fprintf(stderr,"Unable to open input map, or not a .ppm file\n");
  exit(0);
 }

 // Allocate memory for the temporary frame
 fprintf(stderr,"Allocating temp. frame\n");
 map_b=(unsigned char *)calloc(sx*sy*3,sizeof(unsigned char));
 if (map_b==NULL)
 {
  fprintf(stderr,"Out of memory allocating image data\n");
  free(map);
  exit(0);
 }

 srand48((long)time(NULL));		// Initialize random generator from timer
 fprintf(stderr,"Init robot...\n");
 robot=initRobot(map,sx,sy);
 if (robot==NULL)
 {
  fprintf(stderr,"Unable to initialize robot.\n");
  free(map);
  free(map_b);
  exit(0);
 }
 sonar_measurement(robot,map,sx,sy);	// Initial measurements...

 // Initialize particles at random locations
 fprintf(stderr,"Init particles...\n");
 list=NULL;
 initParticles();

 // Done, set up OpenGL and call particle filter loop
 fprintf(stderr,"Entering main loop...\n");
 Win[0]=800;
 Win[1]=800;
 glutInit(&argc, argv);
 initGlut(argv[0]);
 glutMainLoop();

 // This point is NEVER reached... memory leaks must be resolved by OpenGL main loop
 exit(0);

}

void initParticles(void)
{
 /*
   This function creates and returns a linked list of particles
   initialized with random locations (not over obstacles or walls)
   and random orientations.

   There is a utility function to help you find whether a particle
   is on top of a wall.

   Use the global pointer 'list' to keep trak of the *HEAD* of the
   linked list.

   Probabilities should be uniform for the initial set.
 */

 list=NULL;

 /***************************************************************
 // TO DO: Complete this function to generate an initially random
 //        list of particles.
 ***************************************************************/
 for (int i = 0; i < n_particles; i++) {
  struct particle *temp;
  temp = (particle*)malloc(sizeof(particle));
  do {
   temp->x = rand() % sx;
   temp->y = rand() % sy;
  } while (hit(temp, map, sx, sy));
  temp->theta = rand() % 12;
  //temp->theta = temp->theta * 3;
  temp->prob = 1.00000 / n_particles;
  temp->next = list;
  list = temp;
 }


}
double randn (double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }

  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return (mu + sigma * (double) X1);
}

void computeLikelihood(struct particle *p, struct particle *rob, double noise_sigma)
{
 /*
   This function computes the likelihood of a particle p given the sensor
   readings from robot 'robot'.

   Both particles and robot have a measurements array 'measureD' with 16
   entries corresponding to 16 sonar slices.

   The likelihood of a particle depends on how closely its measureD
   values match the values the robot is 'observing' around it. Of
   course, the measureD values for a particle come from the input
   map directly, and are not produced by a simulated sonar.

   Assume each sonar-slice measurement is independent, and if

   error_i = (p->measureD[i])-(sonar->measureD[i])

   is the error for slice i, the probability of observing such an
   error is given by a Gaussian distribution with sigma=20 (the
   noise standard deviation hardcoded into the robot's sonar).

   You may want to check your numbers are not all going to zero...

   This function updates the likelihood for the particle in
   p->prob
 */

 /****************************************************************5000
 // TO DO: Complete this function to calculate the particle's
 //        likelihood given the robot's measurements
 ****************************************************************/
 double error;
 double sum = 1.0;

 for (int i = 0; i< 16; i++) {
    error = p->measureD[i] - rob->measureD[i];
    //printf("ERROR: %f\n", error); 
    //printf("GaussEval: %f\n", GaussEval(error, 20));
    sum *= GaussEval(error, 20);
 }

p->prob = sum;
//printf("p->prob: %f\n", p->prob);
//exit(1);
 
}


void ParticleFilterLoop(void)
{
 /*
    Main loop of the particle filter
 */

  // OpenGL variables. Do not remove
  unsigned char *tmp;
  GLuint texture;
  static int first_frame=1;
  double max;
  struct particle *p,*pmax;
  char line[1024];
  // Add any local variables you need right below.

  if (!first_frame)
  {
   // Step 1 - Move all particles a given distance forward (this will be in
   //          whatever direction the particle is currently looking).
   //          To avoid at this point the problem of 'navigating' the
   //          environment, any particle whose motion would result
   //          in hitting a wall should be bounced off into a random
   //          direction.
   //          Once the particle has been moved, we call ground_truth(p)
   //          for the particle. This tells us what we would
   //          expect to sense if the robot were actually at the particle's
   //          location.
   //
   //          Don't forget to move the robot the same distance!
   //

   /******************************************************************
   // TO DO: Complete Step 1 and test it
   //        You should see a moving robot and sonar figure with
   //        a set of moving particles.
   ******************************************************************/

   move(robot, 1);
   while (hit(robot, map, sx, sy)) {
    move(robot, -2);
    robot->theta += 6;
    move(robot, 2);
   }

   particle *temp1 = list;

   while(temp1!=NULL) {
    move(temp1, 1);
    ground_truth(temp1, map, sx, sy);
    while (hit(temp1, map, sx, sy)) {
     move(temp1, -2);
     temp1->theta += 6;
     move(temp1, 2);
    }
    temp1 = temp1->next;
   }


   // Step 2 - The robot makes a measurement - use the sonar
   sonar_measurement(robot,map,sx,sy);

   // Step 3 - Compute the likelihood for particles based on the sensor
   //          measurement. See 'computeLikelihood()' and call it for
   //          each particle. Once you have a likelihood for every
   //          particle, turn it into a probability by ensuring that
   //          the sum of the likelihoods for all particles is 1.
   struct particle *p;
   p = list;
   struct particle *pp;
   double totalProb;
   double value = 0.2;
   int i = 0;
   while (p != NULL) {	
      computeLikelihood(p, robot, 20.00);
      //printf("p: %i and %f and %f and %f\n", i,p->x,p->y, p->prob);
      totalProb += p->prob;
      p = p->next;
      i++;
   }
  
   p = list;
   //printf("TTToalp: %f\n", totalProb);
   while (p != NULL) {
       //printf("p:  and %f and %f and %f\n",p->x,p->y, p->prob);
       p->prob /= totalProb;
       printf("Ap:  and %f and %f and %f\n",p->x,p->y, p->prob);
       p = p->next;
   }
 
   totalProb = 0.00;

   
  
   /*******************************************************************
   // TO DO: Complete Step 3 and test it
   //        You should see the brightness of particles change
   //        depending on how well they agree with the robot's
   //        sonar measurements. If all goes well, particles
   //        that agree with the robot's position/direction
   //        should be brightest.
   *******************************************************************/

   // Step 4 - Resample particle set based on the probabilities. The goal
   //          of this is to obtain a particle set that better reflect our
   //          current belief on the location and direction of motion
   //          for the robot. Particles that have higher probability will
   //          be selected more often than those with lower probability.
   //
   //          To do this: Create a separate (new) list of particles,
   //                      for each of 'n_particles' new particles,
   //                      randomly choose a particle from  the current
   //                      set with probability given by the particle
   //                      probabilities computed in Step 3.
   //                      Initialize the new particles (x,y,theta)
   //                      from the selected particle.
   //                      Note that particles in the current set that
   //                      have high probability may end up being
   //                      copied multiple times.
   //
   //                      Once you have a new list of particles, replace
   //                      the current list with the new one. Be sure
   //                      to release the memory for the current list
   //                      before you lose that pointer!
   //
   //printf("Ap:  and %f and %f and %f\n",list->x,list->y, list->prob);

    int n = n_particles;
struct particle *head;
double probList[4][n];
double wn[n];
double wr[n];
double sum=0;
double wrsum=0;
head = list;
//setup an array with probability of each particles
//not necessary, but since I am new to C, this works the best for me
for (int i=0; i < n; i++){
probList[0][i] = list->x;
probList[1][i] = list->y;
probList[2][i] = list->theta;
probList[3][i]= list->prob;
//printf("%f\n", list->prob);
sum = sum + list->prob;
list = list->next;
}
//printf("==================================\n");
//normalize the probability
for (int i=0; i < n; i++){
probList[3][i] = n*probList[3][i]/sum;
wn[i] = floor(probList[3][i]);
wr[i] = probList[3][i] - wn[i];
wrsum = wrsum + wr[i];
}
for (int i=0; i < n; i++){
wr[i] = wr[i]/wrsum;
}
//resample and create an array to store the index of the choosen particle
int k=0;
int newProbList[n];
for (int i=0; i<n;i++ ){
//printf("%i\n", i);
for (int j=0; j<wn[i]; j++){
newProbList[k] = i;
//printf("%i\n", i);
//printf("%i\n", newProbList[k]);
k = k+1;
}
}
//printf("================end of i==================\n");
double cs[n];
cs[0] = wr[0];
for (int i=1; i<n; i++){
cs[i] = cs[i-1] + wr[i];
//printf("%f | + %f\n", cs[i],wr[i] );
}
for (int j=k; j<n; j++){
int RandIndex = rand() % n;
while(cs[RandIndex] == 0){
RandIndex = rand() % n;
}
//printf("%i\n", RandIndex);
newProbList[j] = RandIndex;
//printf("%i\n", RandIndex);
//printf("successful\n");
//printf("=============================\n");
//printf("%i\n", j);
}
//printf("===========end of randomness==================\n");
for (int i=0; i<n; i++){
//printf("%i\n", newProbList[0]);
}
//create a new linked list with randomly choosen particles
//head=list
for (int i=0; i<n;i++){
struct particle *nlist;
nlist = (struct particle*)malloc(sizeof(struct particle));
nlist->x = rand() %sx; //probList[0][newProbList[i]];
nlist->y = rand() %sy;// probList[1][newProbList[i]];
nlist->theta = rand() % 12;//probList[2][newProbList[i]];
nlist->prob = probList[3][newProbList[i]];
//exit(0);
nlist->next= list;
list=nlist;
//printf("successful\n");
}
//list = (struct particle*)malloc(sizeof(struct particle));
//printf("===========end of new list==================\n");
//exit(1);
list= head;

   

    
   
   /*******************************************************************
   // TO DO: Complete and test Step 4
   //        You should see most particles disappear except for
   //        those that agree with the robot's measurements.
   //        Hopefully the largest cluster will be on and around
   //        the robot's actual location/direction.
   *******************************************************************/

  }  // End if (!first_frame)

  /***************************************************
   OpenGL stuff
   You DO NOT need to read code below here. It only
   takes care of updating the screen.
  ***************************************************/
  if (RESETflag)	// If user pressed r, reset particles
  {
   deleteList(list);
   list=NULL;
   initParticles();
   RESETflag=0;
  }
  renderFrame(map,map_b,sx,sy,robot,list);

  // Clear the screen and depth buffers
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glEnable(GL_TEXTURE_2D);
  glDisable(GL_LIGHTING);

  glGenTextures( 1, &texture);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glBindTexture( GL_TEXTURE_2D, texture);

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

  glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);

  glTexImage2D (GL_TEXTURE_2D, 0, GL_RGB, sx, sy, 0, GL_RGB, GL_UNSIGNED_BYTE, map_b);

  // Draw box bounding the viewing area
  glBegin (GL_QUADS);
  glTexCoord2f (0.0, 0.0);
  glVertex3f (0.0, 100.0, 0.0);
  glTexCoord2f (1.0, 0.0);
  glVertex3f (800.0, 100.0, 0.0);
  glTexCoord2f (1.0, 1.0);
  glVertex3f (800.0, 700.0, 0.0);
  glTexCoord2f (0.0, 1.0);
  glVertex3f (0.0, 700.0, 0.0);
  glEnd ();

  p=list;
  max=0;
  while (p!=NULL)
  {
   if (p->prob>max)
   {
    max=p->prob;
    pmax=p;
   }
   p=p->next;
  }

  if (!first_frame)
  {
   sprintf(&line[0],"X=%3.2f, Y=%3.2f, th=%3.2f, EstX=%3.2f, EstY=%3.2f, Est_th=%3.2f, Error=%f",robot->x,robot->y,robot->theta,\
           pmax->x,pmax->y,pmax->theta,sqrt(((robot->x-pmax->x)*(robot->x-pmax->x))+((robot->y-pmax->y)*(robot->y-pmax->y))));
   glColor3f(1.0,1.0,1.0);
   glRasterPos2i(5,22);
   for (int i=0; i<strlen(&line[0]); i++)
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,(int)line[i]);
  }

  // Make sure all OpenGL commands are executed
  glFlush();
  // Swap buffers to enable smooth animation
  glutSwapBuffers();

  glDeleteTextures( 1, &texture );

  // Tell glut window to update ls itself
  glutSetWindow(windowID);
  glutPostRedisplay();

  if (first_frame)
  {
   fprintf(stderr,"All set! press enter to start\n");
   gets(&line[0]);
   first_frame=0;
  }
}

/*********************************************************************
 OpenGL and display stuff follows, you do not need to read code
 below this line.
*********************************************************************/
// Initialize glut and create a window with the specified caption
void initGlut(char* winName)
{
    // Set video mode: double-buffered, color, depth-buffered
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

    // Create window
    glutInitWindowPosition (0, 0);
    glutInitWindowSize(Win[0],Win[1]);
    windowID = glutCreateWindow(winName);

    // Setup callback functions to handle window-related events.
    // In particular, OpenGL has to be informed of which functions
    // to call when the image needs to be refreshed, and when the
    // image window is being resized.
    glutReshapeFunc(WindowReshape);   // Call WindowReshape whenever window resized
    glutDisplayFunc(ParticleFilterLoop);   // Main display function is also the main loop
    glutKeyboardFunc(kbHandler);
}

void kbHandler(unsigned char key, int x, int y)
{
 if (key=='r') {RESETflag=1;}
 if (key=='q') {deleteList(list); free(map); free(map_b); exit(0);}
}

void WindowReshape(int w, int h)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();			// Initialize with identity matrix
    gluOrtho2D(0, 800, 800, 0);
    glViewport(0,0,w,h);
    Win[0] = w;
    Win[1] = h;
}
