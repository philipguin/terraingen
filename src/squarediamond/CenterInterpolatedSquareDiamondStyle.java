package squarediamond;
import java.util.Random;



public final class CenterInterpolatedSquareDiamondStyle extends InterpolatedSquareDiamondStyle
{
	public CenterInterpolatedSquareDiamondStyle(Random random)
	{
		super(random);
	}

	@Override
	protected final float nextWeight()
	{
		return .5f;
	}
}
