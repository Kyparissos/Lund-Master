model FurutaPendulum "Furuta pendulum"  
    
    parameter Real pendulumA_start = -Modelica.Constants.pi;   // "-Modelica.Constants.pi"
    
    model ControllerLQR
      Modelica.Blocks.Interfaces.RealInput phi, dphi, theta1, dtheta1, theta2, dtheta2;
      Modelica.Blocks.Interfaces.RealOutput u(start=0);
      parameter Real L[6] = {0, 0, 0, 0, 0, 0};//{100.0, 9.65721819, 103.07438354, -1.40000991, 87.40106372, -3.88918398};
      Real x[6];
    equation
      x = {phi+3.14/2, dphi, theta1+3.14, dtheta1, theta2+3.14, dtheta2};
      u = - L * x;
    end ControllerLQR;
    
    ControllerLQR Controller;
    
    inner Modelica.Mechanics.MultiBody.World world(
      axisLength(displayUnit="mm"),
      axisDiameter(displayUnit="mm"),
      nominalLength(displayUnit="mm") = 0.1)
                                annotation (Placement(transformation(extent={{-60,-60},
              {-40,-40}},   rotation=0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute rotor(
      a(fixed=false),
      w(fixed=true),
      cylinderLength(displayUnit="mm") = 0.015,
      cylinderColor={0,0,0},
      cylinderDiameter(displayUnit="mm") = 0.0605,
      useAxisFlange=true,
      n={0,1,0},
      phi(fixed=true, start=-1.5707963267949)) annotation (Placement(
          transformation(
          origin={0,32},
          extent={{-10,-10},{10,10}},
          rotation=90)));
    Modelica.Mechanics.MultiBody.Joints.Revolute pendulumAxis(
      a(fixed=false),
      w(fixed=true),
      cylinderLength(displayUnit="mm") = 0.005,
      cylinderDiameter(displayUnit="mm") = 0.005,
      cylinderColor={200,200,200},
      useAxisFlange=true,
      n={-1,0,0},
      phi(
        fixed=true,
        start=pendulumA_start,
        displayUnit="rad")) annotation (Placement(transformation(extent={{38,56},
              {58,76}}, rotation=0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder base(
      r(displayUnit="mm") = {0,0.1,0},
      r_shape(displayUnit="mm") = {0,0,0},
      diameter(displayUnit="mm") = 0.06,
      color={155,155,155},
      r_0(displayUnit="mm", fixed=true)) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={0,-26})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder pendulumArm(
      r_shape(displayUnit="mm") = {0,0,0},
      diameter(displayUnit="mm") = 0.005,
      color={200,200,200},
      r(displayUnit="mm") = {0,0.075,0},
      density=3700) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={70,2})));
    Modelica.Mechanics.MultiBody.Parts.Fixed fixed(
      length=0,
      width=0,
      height=0,
      r(displayUnit="mm") = {0,-0.025,-0.1})
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=90,
          origin={0,-52})));
    Modelica.Mechanics.Rotational.Components.Damper pendulumDamper(d=0.000006)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={48,86})));
    Modelica.Mechanics.Rotational.Components.Damper rotorDamper(d=0.009) annotation (
        Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={-24,42})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder pendulumAttachment(
      r_shape(displayUnit="mm") = {0,0,0},
      diameter(displayUnit="mm") = 0.005,
      color={155,155,155},
      r(displayUnit="mm") = {0.043,0,0},
      density=3700) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={14,66})));

    Modelica.Mechanics.Rotational.Sensors.AngleSensor pendulumA
      annotation (Placement(transformation(extent={{78,64},{92,78}})));
    Modelica.Mechanics.Rotational.Sensors.SpeedSensor pendulumW
      annotation (Placement(transformation(extent={{78,80},{92,94}})));
    Modelica.Mechanics.Rotational.Sensors.AngleSensor rotorA
      annotation (Placement(transformation(extent={{-20,-12},{-6,2}})));
    Modelica.Mechanics.Rotational.Sensors.SpeedSensor rotorW
      annotation (Placement(transformation(extent={{-19.937923177416167,2.0},{-6.0620768225838315,16.0}},rotation = 0.0,origin = {0.0,0.0})));
    Modelica.Blocks.Sources.Pulse pulse[3](
      startTime={1,1,1},
      width={1,1,1},
      period={10,10,10},
      amplitude={0.025,0,0})   annotation (Placement(transformation(
          extent={{-6,-6},{6,6}},
          rotation=0,
          origin={32,-32})));
    Modelica.Mechanics.MultiBody.Forces.WorldForce disturbance
      annotation (Placement(transformation(extent={{50,-40},{66,-24}})));
    .Modelica.Mechanics.MultiBody.Parts.BodyCylinder pendulumAttachment2(density = 3700,r(displayUnit = "mm") = {-0.043,0,0},color = {155,155,155},diameter(displayUnit = "mm") = 0.005,r_shape(displayUnit = "mm") = {0,0,0}) annotation(Placement(transformation(extent = {{-107.35765691395954,60.55454850210597},{-131.7407086833764,83.22531007761623}},rotation = 0.0,origin = {0.0,0.0})));
    .Modelica.Mechanics.Rotational.Components.Damper pendulumDamper2(d = 0.000003) annotation(Placement(transformation(extent = {{12.191525884708426,-11.335380787755128},{-12.191525884708426,11.335380787755128}},rotation = -180.0,origin = {-161.00037080667659,94.56069086537136})));
    .Modelica.Mechanics.MultiBody.Joints.Revolute pendulumAxis2(phi(fixed = true,start = pendulumA_start,displayUnit = "rad"),n = {1,0,0},useAxisFlange = true,cylinderColor = {200,200,200},cylinderDiameter(displayUnit = "mm") = 0.005,cylinderLength(displayUnit = "mm") = 0.005,w(fixed = true),a(fixed = false)) annotation(Placement(transformation(extent = {{-148.80884492196816,60.55454850210597},{-173.191896691385,83.22531007761623}},rotation = 0.0,origin = {0.0,0.0})));
    .Modelica.Mechanics.Rotational.Sensors.SpeedSensor pendulumW2 annotation(Placement(transformation(extent = {{-197.57494846080192,87.75946239271829},{-214.64308469939368,103.62899549557548}},rotation = 0.0,origin = {0.0,0.0})));
    .Modelica.Mechanics.Rotational.Sensors.AngleSensor pendulumA2 annotation(Placement(transformation(extent = {{-197.57494846080192,69.62285313231007},{-214.64308469939368,85.49238623516726}},rotation = 0.0,origin = {0.0,0.0})));
    .Modelica.Mechanics.MultiBody.Parts.BodyCylinder pendulumArm2(density = 3700,r(displayUnit = "mm") = {0,0.030,0},color = {200,200,200},diameter(displayUnit = "mm") = 0.005,r_shape(displayUnit = "mm") = {0,0,0}) annotation(Placement(transformation(extent = {{11.335380787755128,-12.191525884708422},{-11.335380787755128,12.191525884708422}},rotation = 90.0,origin = {-187.82172775303513,-0.656507751771727})));
    .Modelica.Mechanics.MultiBody.Forces.WorldForce disturbance2 annotation(Placement(transformation(extent = {{-163.43867598361828,-48.265107060343254},{-182.94511739915174,-30.128497799935047}},rotation = 0.0,origin = {0.0,0.0})));
    .Modelica.Mechanics.Rotational.Sources.Torque torque annotation(Placement(transformation(extent = {{-72.0,22.0},{-52.0,42.0}},origin = {0.0,0.0},rotation = 0.0)));
    equation
    connect(pendulumAxis.frame_b, pendulumArm.frame_a) annotation (Line(
        points={{58,66},{70,66},{70,12}},
        color={95,95,95},
        thickness=0.5,
        smooth=Smooth.None));
    connect(base.frame_b, rotor.frame_a) annotation (Line(
        points={{0,-16},{0,0},{-6.66134e-16,0},{-6.66134e-16,22}},
        color={95,95,95},
        thickness=0.5,
        smooth=Smooth.None));
    connect(base.frame_a, fixed.frame_b) annotation (Line(
        points={{0,-36},{0,-42},{4.44089e-16,-42}},
        color={95,95,95},
        thickness=0.5,
        smooth=Smooth.None));
    connect(pendulumDamper.flange_b, pendulumAxis.support) annotation (Line(
        points={{38,86},{30,86},{30,76},{42,76}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(pendulumDamper.flange_a, pendulumAxis.axis) annotation (Line(
        points={{58,86},{66,86},{66,76},{48,76}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(rotorDamper.flange_a, rotor.axis) annotation (Line(
        points={{-14,42},{-10,42},{-10,32}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(rotorDamper.flange_b, rotor.support) annotation (Line(
        points={{-34,42},{-40,42},{-40,26},{-10,26}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(pendulumAttachment.frame_b, pendulumAxis.frame_a) annotation (Line(
        points={{24,66},{38,66}},
        color={95,95,95},
        thickness=0.5,
        smooth=Smooth.None));
    connect(pendulumAttachment.frame_a, rotor.frame_b) annotation (Line(
        points={{4,66},{0,66},{0,42},{4.44089e-16,42}},
        color={95,95,95},
        thickness=0.5,
        smooth=Smooth.None));
    connect(pendulumW.flange, pendulumAxis.axis) annotation (Line(
        points={{78,87},{74,87},{74,76},{48,76}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(pendulumA.flange, pendulumAxis.axis) annotation (Line(
        points={{78,71},{74,71},{74,76},{48,76}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(rotorW.flange, rotor.axis) annotation (Line(
        points={{-19.937923177416167,9},{-19.937923177416167,32},{-10,32}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(rotorA.flange, rotor.axis) annotation (Line(
        points={{-20,-5},{-24,-5},{-24,32},{-10,32}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(pulse.y,disturbance. force) annotation (Line(
        points={{38.6,-32},{48.4,-32}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(disturbance.frame_b, pendulumArm.frame_b) annotation (Line(
        points={{66,-32},{70,-32},{70,-8}},
        color={95,95,95},
        thickness=0.5,
        smooth=Smooth.None));
    connect(pendulumDamper2.flange_b,pendulumAxis2.support) annotation(Line(points = {{-148.8088449219682,94.56069086537136},{-139.05562421420143,94.56069086537136},{-139.05562421420143,83.22531007761623},{-153.68545527585155,83.22531007761623}},color = {0,0,0}));
    connect(pendulumDamper2.flange_a,pendulumAxis2.axis) annotation(Line(points = {{-173.191896691385,94.56069086537136},{-182.94511739915177,94.56069086537136},{-182.94511739915177,83.22531007761623},{-161.0003708066766,83.22531007761623}},color = {0,0,0}));
    connect(pendulumAttachment2.frame_b,pendulumAxis2.frame_a) annotation(Line(points = {{-131.7407086833764,71.8899292898611},{-148.8088449219682,71.8899292898611}},color = {95,95,95}));
    connect(pendulumW2.flange,pendulumAxis2.axis) annotation(Line(points = {{-197.57494846080186,95.69422894414689},{-192.6983381069185,95.69422894414689},{-192.6983381069185,83.22531007761623},{-161.0003708066766,83.22531007761623}},color = {0,0,0}));
    connect(pendulumA2.flange,pendulumAxis2.axis) annotation(Line(points = {{-197.57494846080186,77.55761968373866},{-192.6983381069185,77.55761968373866},{-192.6983381069185,83.22531007761623},{-161.0003708066766,83.22531007761623}},color = {0,0,0}));
    connect(pendulumAxis2.frame_b,pendulumArm2.frame_a) annotation(Line(points = {{-173.191896691385,71.8899292898611},{-187.82172775303513,71.8899292898611},{-187.82172775303513,10.6788730359834}},color = {95,95,95}));
    connect(disturbance2.frame_b,pendulumArm2.frame_b) annotation(Line(points = {{-182.94511739915177,-39.196802430139165},{-187.82172775303513,-39.196802430139165},{-187.82172775303513,-11.991888539526856}},color = {95,95,95}));
    connect(rotor.frame_b,pendulumAttachment2.frame_a) annotation(Line(points = {{2.220446049250313e-15,42},{2.220446049250313e-15,71.8899292898611},{-107.35765691395954,71.8899292898611}},color = {95,95,95}));
    
    connect(rotorA.phi,Controller.phi);
    connect(pendulumA.phi,Controller.theta1);
    connect(pendulumA2.phi,Controller.theta2);
    connect(rotorW.w,Controller.dphi);
    connect(pendulumW.w,Controller.dtheta1);
    connect(pendulumW2.w,Controller.dtheta2);
    connect(Controller.u,torque.tau);
    
    connect(torque.flange,rotor.axis);
    
    connect(pulse.y,disturbance2.force) annotation(Line(points = {{38.6,-32},{44.6,-32},{44.6,-44},{20.6,-44},{20.6,-68},{-161.48803184206494,-68},{-161.48803184206494,-39.19680243013915}},color = {0,0,127}));
    //rotor.axis
    annotation (
      versionDate="2014-02-04",
      Commands(file="Furuta.mos" "Simulate Furuta pendulum", file="Animate.mos"
          "Animate Furuta pendulum"),
      experiment(NumberOfIntervals=5000, StopTime=10),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}),     graphics),uses(Modelica(version = "4.0.0")));
end FurutaPendulum;
