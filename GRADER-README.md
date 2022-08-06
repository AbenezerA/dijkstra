# HW10 - Sample Grader Classes 

## Setup
No setup steps are required. All grader and associated helper .class files are in your root directory.

## Dijkstra Grader

You can execute the grader via the below command:
```
java DijkstraGrader
```
#### 1. computeEuclideanDistance Test

This test compares your implementation of computeEuclideanDistance with the solution implementation, via the below method.
We are including the implementation of this test functionality to help you interpret your grader results and debug from there.

```
  private static boolean testComputeEuclideanDistance(){
    Dijkstra studentCode = new Dijkstra();
    DijkstraSolution solutionCode = new DijkstraSolution();

    final double[][] distances = {{0.0, 0.0, 0.0, 0.0},
                            {10.0, 10.0, 10.0, 10.0},
                            {-10.0, -10.0, -10.0, -10.0},
                            {10.0, 20.0, 30.0, 40.0},
                            {-10.0, -20.0, -30.0, -40.0},
                            {10.0, -20.0, -30.0, 40.0}
                           };
    for(int i = 0; i < distances.length; i++){
      double studentDistance = studentCode.computeEuclideanDistance(distances[i][0], distances[i][1], 
                                                                    distances[i][2], distances[i][3]);
      double solutionDistance = solutionCode.computeEuclideanDistance(distances[i][0], distances[i][1], 
                                                                    distances[i][2], distances[i][3]);
      if(Math.abs(solutionDistance - studentDistance) > EPSILON){
        return false;
      }
    }
    return true;
  } 
```
#### 2. computeAllEuclideanDistances Test
This test compares your implementation of computeAllEuclideanDistances with the solution implementation, via the below method.
We are including the implementation of this test functionality to help you interpret your grader results and debug from there.

```
  private static boolean testComputeAllEuclideanDistances() throws IOException{
    Dijkstra studentCode = new Dijkstra();
    DijkstraSolution solutionCode = new DijkstraSolution();

    Map<String, Vertex> studentMap = studentCode.getVertexMap();
    Map<String, Vertex> solutionMap = solutionCode.getVertexMap();
    
    for(String name : solutionMap.keySet()){
      List<Edge> solutionEdges = solutionMap.get(name).adjacentEdges;
      List<Edge> studentEdges = studentMap.get(name).adjacentEdges;
      
      for(int i = 0; i < solutionEdges.size(); i++){
        Edge solutionEdge = solutionEdges.get(i);
        Edge studentEdge = studentEdges.get(i);
        if(!solutionEdge.equals(studentEdge)){
          return false;
        }
      }
    }
    return true;
  }
  
```
#### 3. Dijkstra Test: NewYork to SanFrancisco.
This test compares your full implementation with the solution implementation, via the below method.
We are including the implementation of this test functionality to help you interpret your grader results and debug from there.

```
    Dijkstra studentCode = new Dijkstra();
    DijkstraSolution solutionCode = new DijkstraSolution();

    String startCity = "NewYork"
    String endCity = "SanFrancisco";
    
    try{
      if(testDijkstra(studentCode, solutionCode, startCity, endCity)){
        System.out.println("Test 1 success: NewYork to SanFrancisco. ");
        correctTest.add(0);
      }
      else{
        errorLog.add("-6 Failed Test 1: NewYork to SanFrancisco. ");
        totalScore -= 6.0;
      }
    }
    catch(Throwable e){
      errorLog.add("-6 Test 1 threw an exception. ");
      totalScore -= 6.0;
    }
```

The testDijkstra method that's invoked above essentially does the following to test your paths:

```
    try{
      studentCode.doDijkstra(s);
      solutionCode.doDijkstra(s);
      List<Edge> studentPath = studentCode.getDijkstraPath(s, t);
      List<Edge> solutionPath = solutionCode.getDijkstraPath(s, t);
      if(testPath(studentPath, solutionPath)){
          return true;
      }
      else{
        return false;
      }
    }
    catch(Throwable e){
      throw new Exception();
    }
```
Just a reminder that these tests are not exhaustive (neither of the requirements nor of possible edge cases). During grading, your doDijkstra, getDijkstraPath implementations will be seperately graded, as well as the testing of further start and end locations. As such, please perform additional testing on your classes. These released tests make up about a third of the the test cases you will be graded with. 