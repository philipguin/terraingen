package erosion;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.Random;

import src.Array2D;
import src.IModifier;

public class RiverModifier implements IModifier<Array2D>
{
	private final Random random;
	private final float chanceOfRiverSpawn;
	
	public RiverModifier(Random random, float chanceOfRiverSpawn)
	{
		this.random = random;
		this.chanceOfRiverSpawn = chanceOfRiverSpawn;
	}
	
	private static final class Point
	{
		public final int x, y;
		
		public Point(int x, int y)
		{
			this.x = x;
			this.y = y;
		}
	}
	
	private final Deque<Point> stack = new ArrayDeque<Point>();
	
	@Override
	public void modify(Array2D values)
	{
		for (int r = (int)Math.ceil(chanceOfRiverSpawn * values.width * values.height); r > 0; --r)
		{
			int x, y;
			float currentHeight;
			
			do
			{
				x = random.nextInt(values.width);
				y = random.nextInt(values.height);
				currentHeight = values.get(x, y);
			}
			while (currentHeight < .7f);
			
			int direction;
			
			while (currentHeight >= .45f)
			{
				stack.add(new Point(x, y));
				
				direction = -1;
				
				if (x > 0 && values.get(x - 1, y) < currentHeight)
				{
					direction = 0;
					currentHeight = values.get(x - 1, y);
				}
				
				if (y > 0 && values.get(x, y - 1) < currentHeight)
				{
					direction = 1;
					currentHeight = values.get(x, y - 1);
				}
				
				if (x < values.width - 1 && values.get(x + 1, y) < currentHeight)
				{
					direction = 2;
					currentHeight = values.get(x + 1, y);
				}
				
				if (y < values.height - 1 && values.get(x, y + 1) < currentHeight)
				{
					direction = 3;
					currentHeight = values.get(x, y + 1);
				}
				
				if (direction == -1)
					break;
					
				switch (direction)
				{
				case 0: --x; break;
				case 1: --y; break;
				case 2: ++x; break;
				case 3: ++y; break;
				}
			}
		}
		
		for (Point p : stack)
		{
			values.set(p.x, p.y, Math.max(0f, values.get(p.x, p.y) - .05f));

			if (p.x > 0) values.set(p.x - 1, p.y, Math.max(0f, values.get(p.x - 1, p.y) - .025f));
			if (p.y > 0) values.set(p.x, p.y - 1, Math.max(0f, values.get(p.x, p.y - 1) - .025f));
			if (p.x < values.width  - 1) values.set(p.x + 1, p.y, Math.max(0f, values.get(p.x + 1, p.y) - .025f));
			if (p.y < values.height - 1) values.set(p.x, p.y + 1, Math.max(0f, values.get(p.x, p.y + 1) - .025f));
		}
		
		stack.clear();
	}

}
