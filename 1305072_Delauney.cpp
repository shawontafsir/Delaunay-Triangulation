#include <bits/stdc++.h>
#include <windows.h>
#include <glut.h>
#define PI 3.14159265

using namespace std;

class Point{
public:
    double x;
    double y;
    double angle;

    Point(){}
    Point(double x, double y){this->x = x; this->y = y;}

    bool operator < (const Point& rhs) const
    {
        return (angle < rhs.angle);
    }
    bool operator == (const Point& rhs) const
    {
        return (x==rhs.x && y==rhs.y);
    }
    bool operator > (const Point& p) const
    {
        if(y==p.y) return x>p.x;
        else return y>p.y;
    }
};

class Edge{
public:
    Point bottom, top;
    Edge(){}
    Edge(Point p1, Point p2){
        if(p1.y<p2.y){
            bottom = p1;
            top = p2;
        }
        else if(p1.y>p2.y){
            bottom = p2;
            top = p1;
        }
        else{
            if(p1.x<p2.x){
                bottom = p1;
                top = p2;
            }
            else{
                bottom = p2;
                top = p1;
            }
        }
    }

    bool operator < (const Edge& edge) const
    {
        if(top==edge.top) return edge.bottom>bottom;
        return top>edge.top;
    }
};

class Triangle{
public:
    Point p1;
    Point p2;
    Point p3;
    Point cirCenter;

    Triangle(){}
    Triangle(Point p1, Point p2, Point p3){
        this->p1 = p1;
        this->p2 = p2;
        this->p3 = p3;
    }
};

vector<Point> glpoints;
vector<Triangle> glTau;
int glpointssize, glTausize, voronoiSwitch = 0;

void drawAxes(){
    glColor3f(0,0,0);
	glBegin(GL_LINES);
    {
        glVertex2f(200*.09, 0*.09);
        glVertex2f(-200*.09, 0*.09);
    }
    glEnd();
    glBegin(GL_LINES);
    {
        glVertex2f(0,200*.09);
        glVertex2f(0,-200*.09);
    }
    glEnd();
}

double angle_bet_points(Point a, Point b){
    double t;
    if(a.x!=b.x && a.y!=b.y) t = atan((double)(b.y-a.y)/(double)(b.x-a.x));
    else if(a.x==b.x && a.y==b.y) t = 0;
    else if(a.x==b.x) t = (b.y>a.y?1:3)*PI/2;
    else t = (b.x>a.x?0:1)*PI;

    if(((b.x-a.x)<0 && (b.y-a.y)>0) || ((b.y-a.y)<0 && (b.x-a.x)<0)) t+=PI;
    else if((b.y-a.y)<0 && (b.x-a.x)>0) t+=2*PI;

    return t;
}

double area(Point p1, Point p2, Point p3){
    return abs((p1.x*(p2.y-p3.y) + p2.x*(p3.y-p1.y) + p3.x*(p1.y-p2.y))/2.0);
}

bool isInside(Triangle t, Point p){
    double a = area(t.p1, t.p2, t.p3);
    double a1 = area(t.p1, t.p2, p);
    double a2 = area(p, t.p2, t.p3);
    double a3 = area(t.p1, p, t.p3);
    if(abs(a-a1-a2-a3)<pow(10, -10)) return true;
    else return false;
}

bool hasEdge(Triangle t, Point p1, Point p2){
    if((t.p1==p1 && t.p2==p2) || (t.p1==p2 && t.p2==p1)) return true;
    else if((t.p2==p1 && t.p3==p2) || (t.p2==p2 && t.p3==p1)) return true;
    else if((t.p3==p1 && t.p1==p2) || (t.p3==p2 && t.p1==p1)) return true;
    else return false;
}

bool onLine(Point p1, Point p2, Point p){
    double a = (p.y-p1.y)/(p1.y-p2.y) - (p.x-p1.x)/(p1.x-p2.x);

    if(abs(a)<pow(10, -10)) return true;
    else return false;
}

Point circumCenter1(Point p1, Point p2, Point p3){
    Point center;
    center.x = ((p3.y-p1.y)/2-(p1.x-p2.x)/(p1.y-p2.y)*(p1.x+p2.x)/2+(p2.x-p3.x)/(p2.y-p3.y)*(p2.x+p3.x)/2)/((p2.x-p3.x)/(p2.y-p3.y)-(p1.x-p2.x)/(p1.y-p2.y));
    center.y = ((p3.x-p1.x)/2-(p1.y-p2.y)/(p1.x-p2.x)*(p1.y+p2.y)/2+(p2.y-p3.y)/(p2.x-p3.x)*(p2.y+p3.y)/2)/((p2.y-p3.y)/(p2.x-p3.x)-(p1.y-p2.y)/(p1.x-p2.x));

    return center;
}

Point circumCenter(Point p1, Point p2, Point p3) {
    double xysq1, xysq2, xysq3;
    double cirA, cirB, cirC, cirD;
    Point center;
    xysq1 = p1.x*p1.x+p1.y*p1.y;
    xysq2 = p2.x*p2.x+p2.y*p2.y;
    xysq3 = p3.x*p3.x+p3.y*p3.y;

    cirA = p1.x*(p2.y-p3.y) - p1.y*(p2.x-p3.x) + p2.x*p3.y - p3.x*p2.y;
    cirB = xysq1*(p3.y-p2.y) + xysq2*(p1.y-p3.y) + xysq3*(p2.y-p1.y);
    cirC = xysq1*(p2.x-p3.x) + xysq2*(p3.x-p1.x) + xysq3*(p1.x-p2.x);
    cirD = xysq1*(p3.x*p2.y-p2.x*p3.y) + xysq2*(p1.x*p3.y-p3.x*p1.y) + xysq3*(p2.x*p1.y-p1.x*p2.y);

    center.x = -cirB / (2*cirA);
    center.y = -cirC / (2*cirA);

    return center;
}

bool isIllegal(Point p1, Point p2, Point p3, Point pk){
    Point center = circumCenter(p1, p2, p3);
    long double radius = sqrt(pow(p2.x-center.x, 2)+pow(p2.y-center.y, 2));
    long double dist = sqrt(pow(pk.x-center.x, 2)+pow(pk.y-center.y, 2));
    //cout<<"dist and radius : "<<dist<<"  "<<radius<<endl;
    if(radius>dist) return true;
    else return false;
}

void legalizeEdge(Point pr, Point p1, Point p2, vector<Triangle> &tau){
    vector<unsigned int> arr;
    Point pk(1000,100000);
    for(unsigned int j=0;j<tau.size();j++){
        if(hasEdge(tau[j], p1, p2)){
            if(pr == tau[j].p1 || pr == tau[j].p2 || pr == tau[j].p3) arr.insert(arr.begin(), j);
            else{
                arr.push_back(j);
                if((p1 == tau[j].p1 && p2 == tau[j].p2) || (p2 == tau[j].p1 && p1 == tau[j].p2))  pk=tau[j].p3;
                else if((p1 == tau[j].p2 && p2 == tau[j].p3) || (p2 == tau[j].p2 && p1 == tau[j].p3)) pk=tau[j].p1;
                else pk=tau[j].p2;
            }
        }

        if(arr.size()==2) break;
    }

    if(isIllegal(p1, pr, p2, pk)){
        tau.push_back(Triangle(pr, p1, pk));
        tau.push_back(Triangle(pr, p2, pk));

        if(arr[0]>arr[1]) {
            tau.erase(tau.begin()+arr[0]);
            tau.erase(tau.begin()+arr[1]);
        }
        else{
            tau.erase(tau.begin()+arr[1]);
            tau.erase(tau.begin()+arr[0]);
        }

        legalizeEdge(pr, p1, pk, tau);
        legalizeEdge(pr, pk, p2, tau);
    }
}

vector<Triangle> delaunayTriangulation(vector<Point> points, Point p0, int N){
    vector<Triangle> tau;

    for(int i=0;i<N;i++) points[i].angle = angle_bet_points(p0, points[i]);
    sort(points.begin(),points.end());
    Point p1(points[1].x-100*cos(points[1].angle-.1-PI), points[1].y-100*sin(points[1].angle-.1-PI));
    Point p2(points[N-1].x+100*cos(2*PI-points[N-1].angle-.1), points[N-1].y-100*sin(2*PI-points[N-1].angle-.1));
    tau.push_back(Triangle(p0, p1, p2));

    random_shuffle(points.begin()+1, points.end());
    for(unsigned int j=0;j<points.size();j++) cout<<points[j].x<<" "<<points[j].y<<endl;

    for(int i=1;i<N;i++){
        for(unsigned int j=0;j<tau.size();j++){
            if(isInside(tau[j], points[i])){
                if(onLine(tau[j].p1, tau[j].p2, points[i])){
                    Triangle temp = tau[j];
                    //cout<<"online1"<<endl;
                    Point pk(1000, 1000);
                    unsigned int ind;
                    for(unsigned int k=0;k<tau.size();k++){
                        if(k==j) continue;
                        if(hasEdge(tau[k], tau[j].p1, tau[j].p2)){
                            ind = k;
                            if((temp.p1 == tau[k].p1 && temp.p2 == tau[k].p2)||(temp.p2 == tau[k].p1 && temp.p1 == tau[k].p2))  pk=tau[j].p3;
                            else if((temp.p1 == tau[k].p2 && temp.p2 == tau[k].p3) || (temp.p2 == tau[k].p2 && temp.p1 == tau[k].p3)) pk=tau[j].p1;
                            else pk=tau[j].p2;

                            break;
                        }
                    }

                    tau.push_back(Triangle(temp.p1, temp.p3, points[i]));
                    tau.push_back(Triangle(temp.p2, temp.p3, points[i]));
                    tau.push_back(Triangle(temp.p1, pk, points[i]));
                    tau.push_back(Triangle(temp.p2, pk, points[i]));

                    if(ind>j) {
                        tau.erase(tau.begin()+ind);
                        tau.erase(tau.begin()+j);
                    }
                    else{
                        tau.erase(tau.begin()+j);
                        tau.erase(tau.begin()+ind);
                    }

                    cout<<"1";
                    legalizeEdge(points[i], temp.p1, pk, tau);cout<<"2";
                    legalizeEdge(points[i], pk, temp.p2, tau);cout<<"3";
                    legalizeEdge(points[i], temp.p2, temp.p3, tau);cout<<"4";
                    legalizeEdge(points[i], temp.p3, temp.p1, tau);cout<<"5";
                }
                else if(onLine(tau[j].p2, tau[j].p3, points[i])){
                    Triangle temp = tau[j];
                    //cout<<"online2"<<endl;
                    Point pk(1000, 1000);
                    unsigned int ind;
                    for(unsigned int k=0;k<tau.size();k++){
                        if(k==j) continue;
                        if(hasEdge(tau[k], temp.p2, temp.p3)){
                            ind = k;
                            if((temp.p2 == tau[k].p1 && temp.p3 == tau[k].p2)||(temp.p3 == tau[k].p1 && temp.p2 == tau[k].p2)) pk=tau[k].p3;
                            else if((temp.p2 == tau[k].p2 && temp.p3 == tau[k].p3)||(temp.p3 == tau[k].p2 && temp.p2 == tau[k].p3)) pk=tau[k].p1;
                            else pk=tau[k].p2;

                            break;
                        }
                    }

                    tau.push_back(Triangle(temp.p2, temp.p1, points[i]));
                    tau.push_back(Triangle(temp.p3, temp.p1, points[i]));
                    tau.push_back(Triangle(temp.p2, pk, points[i]));
                    tau.push_back(Triangle(temp.p3, pk, points[i]));

                    if(ind>j) {
                        tau.erase(tau.begin()+ind);
                        tau.erase(tau.begin()+j);
                    }
                    else{
                        tau.erase(tau.begin()+j);
                        tau.erase(tau.begin()+ind);
                    }

                    cout<<"1";
                    legalizeEdge(points[i], temp.p2, pk, tau);cout<<"2";
                    legalizeEdge(points[i], pk, temp.p3, tau);cout<<"3";
                    legalizeEdge(points[i], temp.p3, temp.p1, tau);cout<<"4";
                    legalizeEdge(points[i], temp.p1, temp.p2, tau);cout<<"5";
                }
                else if(onLine(tau[j].p3, tau[j].p1, points[i])){
                    Triangle temp = tau[j];
                    //cout<<"tau index : "<<j<<" point : "<<points[i].x<<"  "<<points[i].y<<endl;
                    cout<<tau[j].p1.x<<" "<<tau[j].p1.y<<"   ";
                    cout<<tau[j].p2.x<<" "<<tau[j].p2.y<<"   ";
                    cout<<tau[j].p3.x<<" "<<tau[j].p3.y<<endl;
                    //cout<<"online3"<<endl;
                    Point pk(1000, 1000);
                    unsigned int ind;
                    for(unsigned int k=0;k<tau.size();k++){
                        if(k==j) continue;
                        if(hasEdge(tau[k], temp.p3, temp.p1)){
                            ind = k;
                            if((temp.p3 == tau[k].p1 && temp.p1 == tau[k].p2)||(temp.p1 == tau[k].p1 && temp.p3 == tau[k].p2)) pk=tau[k].p3;
                            else if((temp.p3 == tau[k].p2 && temp.p1 == tau[k].p3)||(temp.p1 == tau[k].p2 && temp.p3 == tau[k].p3)) pk=tau[k].p1;
                            else pk=tau[k].p2;

                            break;
                        }
                    }

                    cout<<tau[ind].p1.x<<" "<<tau[ind].p1.y<<"   ";
                    cout<<tau[ind].p2.x<<" "<<tau[ind].p2.y<<"   ";
                    cout<<tau[ind].p3.x<<" "<<tau[ind].p3.y<<"   pk : "<<pk.x<<" "<<pk.y<<"  "<<ind<<" "<<j<<endl;

                    if(ind>j) {
                        tau.erase(tau.begin()+ind);
                        tau.erase(tau.begin()+j);
                    }
                    else{
                        tau.erase(tau.begin()+j);
                        tau.erase(tau.begin()+ind);
                    }

                    tau.push_back(Triangle(temp.p3, temp.p2, points[i]));
                    tau.push_back(Triangle(temp.p1, temp.p2, points[i]));
                    tau.push_back(Triangle(temp.p1, pk, points[i]));
                    tau.push_back(Triangle(temp.p3, pk, points[i]));

                    cout<<"1";
                    legalizeEdge(points[i], temp.p3, pk, tau);
                    legalizeEdge(points[i], pk, temp.p1, tau);
                    legalizeEdge(points[i], temp.p1, temp.p2, tau);
                    legalizeEdge(points[i], temp.p2, temp.p3, tau);
                }
                else{
                    //cout<<"inside "<<points[i].x<<" "<<points[i].y<<endl;
                    Triangle tem = tau[j];
                    tau.push_back(Triangle(tem.p1, tem.p2, points[i]));
                    tau.push_back(Triangle(tem.p2, tem.p3, points[i]));
                    tau.push_back(Triangle(tem.p3, tem.p1, points[i]));
                    tau.erase(tau.begin()+j);

                    legalizeEdge(points[i], tem.p1, tem.p2, tau);
                    legalizeEdge(points[i], tem.p2, tem.p3, tau);
                    legalizeEdge(points[i], tem.p3, tem.p1, tau);
                }

                break;
            }
        }
    }

    for(unsigned int i=0;i<tau.size();i++)
        if(isInside(tau[i], p1) || isInside(tau[i], p2)) {tau.erase(tau.begin()+i);i--;}

    return tau;
}

void voronoi(vector<Triangle> tau){
    glColor3f(0,0,1);
    for(unsigned int i=0;i<tau.size();i++){
        tau[i].cirCenter = circumCenter(tau[i].p1, tau[i].p2, tau[i].p3);
        glBegin(GL_POINTS);
        {
            glVertex2f(tau[i].cirCenter.x*.09, tau[i].cirCenter.y*.09);
        }
        glEnd();
    }

    map<Edge, vector<Triangle> > adj;

    for(unsigned int i=0;i<tau.size();i++){
        Edge e1(tau[i].p1, tau[i].p2);
        adj[e1].push_back(tau[i]);
        Edge e2(tau[i].p2, tau[i].p3);
        adj[e2].push_back(tau[i]);
        Edge e3(tau[i].p3, tau[i].p1);
        adj[e3].push_back(tau[i]);
    }

    map<Edge, vector<Triangle> >::iterator itr;

    for(itr=adj.begin();itr!=adj.end();itr++){
        Edge e = (*itr).first;
        if(adj[e].size()>1){
            glBegin(GL_LINES);
            {
                glVertex2f(adj[e][0].cirCenter.x*.09, adj[e][0].cirCenter.y*.09);
                glVertex2f(adj[e][1].cirCenter.x*.09, adj[e][1].cirCenter.y*.09);
            }
            glEnd();
        }
        else{
            Point p1,p2;
            p1.x = adj[e][0].cirCenter.x;
            p1.y = adj[e][0].cirCenter.y;
            if(!isInside(adj[e][0], p1)) continue;

            double a = (e.top.x+e.bottom.x)/2, b = (e.top.y+e.bottom.y)/2;

            p2.x = 2*a-p1.x;
            p2.y = 2*b-p1.y;

            glBegin(GL_LINES);
            {
                glVertex2f(p1.x*.09, p1.y*.09);
                glVertex2f(p2.x*.09, p2.y*.09);
            }
            glEnd();
        }
    }
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){
		case '1':
			voronoiSwitch = 1-voronoiSwitch;
			glutPostRedisplay();
			break;
		default:
			break;
	}
}

void display(){
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(1,1,1,0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);

	drawAxes();

	glPointSize(5);

	for(int i=0;i<glpointssize;i++){
        glColor3f(1,0,0);

        glBegin(GL_POINTS);
        {
            glVertex2f(glpoints[i].x*.09, glpoints[i].y*.09);
        }
        glEnd();
    }

    for(int i=0;i<glTausize;i++){
        glBegin(GL_LINES);
        {
            glVertex2f(glTau[i].p1.x*.09, glTau[i].p1.y*.09);
            glVertex2f(glTau[i].p2.x*.09, glTau[i].p2.y*.09);
        }
        glEnd();

        glBegin(GL_LINES);
        {
            glVertex2f(glTau[i].p2.x*.09, glTau[i].p2.y*.09);
            glVertex2f(glTau[i].p3.x*.09, glTau[i].p3.y*.09);
        }
        glEnd();

        glBegin(GL_LINES);
        {
            glVertex2f(glTau[i].p1.x*.09, glTau[i].p1.y*.09);
            glVertex2f(glTau[i].p3.x*.09, glTau[i].p3.y*.09);
        }
        glEnd();
    }

    if(voronoiSwitch==1) voronoi(glTau);

	glutSwapBuffers();
}

int main(int argc, char **argv)
{
    Point p0(-90,-90);
    vector<Point> points;
    ifstream in("C:\\Users\\USER\\Favorites\\Desktop\\1305072_Delauney\\1305072_input1.txt");
    if (in.is_open())
    {
        string line;
        getline(in,line);
        int N=atoi(line.c_str());

        for(int i=0;i<N;i++)
        {
            getline(in,line);
            stringstream ss(line);
            string s;
            Point p;
            ss>>s;
            p.x = atof(s.c_str());
            ss>>s;
            p.y = atof(s.c_str());

            if(p.y>p0.y) p0 = p;
            else if(p.y == p0.y && p.x>p0.x) p0 = p;

            points.push_back(p);
        }
        in.close();
    }

    vector<Triangle> tau = delaunayTriangulation(points, p0, points.size());

    glpoints = points;
    glpointssize = points.size();
    glTau = tau;
    glTausize = tau.size();

    glutInit(&argc,argv);
	glutInitWindowSize(600, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
	glutCreateWindow("My OpenGL Program");
	glEnable(GL_DEPTH_TEST);

	glutDisplayFunc(display);
	glutKeyboardFunc(keyboardListener);
	glutMainLoop();

    return 0;
}
