package squarediamond;
import java.util.Random;



public final class UniformRandomlyInterpolatedSquareDiamondStyle extends InterpolatedSquareDiamondStyle
{
	public UniformRandomlyInterpolatedSquareDiamondStyle(Random random)
	{
		super(random);
	}

	@Override
	protected final float nextWeight()
	{
		return random.nextFloat();
	}
}
