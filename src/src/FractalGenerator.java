package src;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.util.Arrays;
import java.util.Random;

import squarediamond.SquareDiamondArray2DPopulator;
import squarediamond.UniformBiasedRandomStyle;
import squarediamond.UniformRandomlyInterpolatedSquareDiamondStyle;
import tectonics.Lithosphere;
import colormap.CachedMapper;
import colormap.Mapper;
import colormap.MountainousColorMap;
import draw.IDrawable;
import erosion.DropletErosionStyle;
import erosion.ErosionModifier;
import frameworks.DrawStepper;
import frameworks.IStepper;
import frameworks.StepLooper;


public class FractalGenerator extends KeyAdapter implements IStepper, IDrawable<Graphics2D>
{
	private final Array2D baseValues;
	private final Mapper<Float, Color> colorMap;
	private final IPopulator<? super Array2D> populator;
	private final IModifier<Array2D> modifier;
	private final int scale;
	
	private Array2D derivedValues, toDraw;
	
	public FractalGenerator(Array2D values, Mapper<Float, Color> colorMap, IPopulator<? super Array2D> populator, IModifier<Array2D> modifier, int scale)
	{
		this.baseValues = values;
		this.colorMap = colorMap;
		this.populator = populator;
		this.modifier = modifier;
		this.scale = scale;
	}
	
	private volatile boolean shouldPopulate = true, shouldModify = false, shouldSwap = true;
	
	@Override
	public void keyTyped(KeyEvent e)
	{
		char c = e.getKeyChar();
		
		if (c == ' ')
			shouldPopulate = true;
		else if (c == 'm')
			shouldModify = true;
		else if (c == 's')
			shouldSwap = true;
	}
	
	@Override
	public void step()
	{
		if (shouldPopulate)
		{
			populator.populate(baseValues);
			shouldPopulate = false;
			//shouldModify = true;
		}
		
		if (shouldModify)
		{
			//derivedValues = modifier.derive(baseValues);
			modifier.modify(baseValues);
			shouldModify = false;
			
			if (toDraw != baseValues)
				toDraw = derivedValues;
		}
		
		if (shouldSwap)
		{
			if (toDraw == baseValues)
				toDraw = derivedValues;
			else
				toDraw = baseValues;
			
			shouldSwap = false;
		}
	}
	
	//private Component component;

	@Override
	public void draw(Graphics2D g)
	{
		int i, j, x, y;
		
		g.setColor(Color.BLACK);
		g.fillRect(0, 0, toDraw.width * scale, toDraw.height * scale);
		
		int halfScale = scale / 2;
		/*Point mousePos = MouseInfo.getPointerInfo().getLocation();
		Point screenPos = component.getLocationOnScreen();
		
		int mouseX = mousePos.x - screenPos.x;
		int mouseY = mousePos.y - screenPos.y;
		int minDistance = 200;*/
		
		for (i = 0, x = halfScale; i < toDraw.width;  ++i, x += scale)
		for (j = 0, y = halfScale; j < toDraw.height; ++j, y += scale)
		{
			g.setColor(colorMap.map(toDraw.get(i, j)));

			/*if (Math.abs(x - mouseX) < minDistance && Math.abs(y - mouseY) < minDistance)
			{
				float distance = (float)Math.hypot(x - mouseX, y - mouseY);
				
				if (distance <= minDistance)
				{
					double angle = Math.atan2(y - mouseY, x - mouseX);
					
					g.fillRect(
							x + (int)(Math.cos(angle) * distance * (1f - (float)Math.pow(distance / minDistance, .2)) + .5d) - halfScale - 1,
							y + (int)(Math.sin(angle) * distance * (1f - (float)Math.pow(distance / minDistance, .2)) + .5d) - halfScale - 1,
							scale + 2,
							scale + 2);
					
					continue;
				}
			}*/
			
			g.fillRect(x - halfScale, y - halfScale, scale, scale);
		}
	}
	
	public static void main(String[] args)
	{
		int size = 8;
		Array2D values = new Array2D((2 << size) + 1, (2 << size) + 1);
		
		Random random = new Random();
		
		Mapper<Float, Color> colorMap = 
			new CachedMapper<Color>(
					new MountainousColorMap(),//new SimpleHueColorMap(2f, 0f, .3f, 1f),
					0x10000);
		
		@SuppressWarnings("unchecked")
		IPopulator<? super Array2D> populator = 
			//new PopulatorModifier<Array2D>(
			new WeightedArray2DAdder(
				Arrays.<IPopulator<Array2D>>asList(
					new SquareDiamondArray2DPopulator(new UniformBiasedRandomStyle(
							random,
							new UniformRandomlyInterpolatedSquareDiamondStyle(
									random/*,
									Arrays.asList(1f, 1f, .8f, .8f, .8f, .8f, .6f, .6f, .4f)*/),
							Arrays.asList(1f,1f, 1f, 1f, .5f, .1f, 0f))
					),
					new SquareDiamondArray2DPopulator(new UniformBiasedRandomStyle(
							random,
							new UniformRandomlyInterpolatedSquareDiamondStyle(
									random/*,
									Arrays.asList(1f, 1f, .8f, .8f, .8f, .8f, .6f, .6f, .4f)*/),
							Arrays.asList(1f,1f, .2f, .2f, .1f, .05f, 0f))
					),
					new PerlinNoiseArray2DPopulator(random, 32, 0f, 1f),
					new PerlinNoiseArray2DPopulator(random, 16, 0f, 1f),
					new PerlinNoiseArray2DPopulator(random, 8,  0f, 1f),
					new PerlinNoiseArray2DPopulator(random, 4,  0f, 1f)
				),
				Arrays.asList(80f, 40f, 8f, 4f, 2f, 1f),
				1f
			);

		@SuppressWarnings("unchecked")
		IModifier<Array2D> modifier = 
			//new RiverModifier(random, .001f)),
			new MultiModifier<Array2D>(
				/*new ErosionModifier(random, 10000, 2f, 0.001f, .9f, 0.05f, .04f, 0.05f, 20f,
					new RainErosionStyle(
							random,
							new WindElevationDerivator(0f, .001f, 0f),
							10f
						)
				),*/
				new ErosionModifier(random, 10000, 2f, 0.001f, .9f, 0.05f, .04f, 0.05f, 15f,
						new DropletErosionStyle(random, 0f, 0f)
				)
			);
		
		@SuppressWarnings("unused")
		IModifier<Array2D> tectonics =  new Lithosphere(
				random,
				5,
				.5f,
				1000, 
				.4f,
				(1 << size) / 7,
				.5f,
				2);
		
		final FractalGenerator gen = new FractalGenerator(
				values,
				colorMap,
				populator,
				//tectonics,
				modifier,
				1 << (9 - size));
		
		final DrawStepper drawStepper = new DrawStepper("Terrain Gen", values.width * gen.scale, values.height * gen.scale, gen);
		drawStepper.getCanvas().addKeyListener(gen);
		//gen.component = window.getFrame();
		
		StepLooper looper = new StepLooper(new IStepper()
		{
			@Override
			public void step()
			{
				gen.step();
				drawStepper.step();
			}
		});
		looper.run();
	}
}
